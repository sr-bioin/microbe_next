#!/usr/bin/env nextflow

/*
 * Bacterial Whole Genome Sequencing Analysis Pipeline
 * Includes QC, Assembly, and Annotation
 */

// Define parameters at the start of the script
params.reads = "data/*_R{1,2}.fastq.gz"
params.outdir = "results"
params.max_memory = '8.GB'
params.max_cpus = 4
params.max_time = '4.h'

// Print pipeline parameters
log.info """
    Bacterial WGS Pipeline
    ======================
    Reads          : ${params.reads}
    Outdir  	   : ${params.outdir}
    
    Max resources:
    CPUs          : ${params.max_cpus}
    Memory        : ${params.max_memory}
    Time          : ${params.max_time}
    """.stripIndent()

// Define channels
reads = Channel.fromFilePairs(params.reads, size: 2)

process fastqc {
    tag "${sample_id}"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path("*.{html,zip}"), emit: reports  // Combined output
    path("*.html"), emit: html           // Separate HTML output
    path("*.zip"), emit: zip             // Separate ZIP output
    
    script:
    """
    fastqc -q -t ${task.cpus} ${reads}
    """
}

process multiqc {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
    input:
    path fastqc_reports
    
    output:
    path("multiqc_report.html"), emit: report
    
    script:
    """
    mkdir -p fastqc_reports
    cp *.html *.zip fastqc_reports/ 2>/dev/null || :
    multiqc fastqc_reports/ -n multiqc_report.html
    """
}

process trimmomatic {
    tag "${sample_id}"
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("*_paired.fq.gz"), emit: trimmed_reads
    path("*.log"), emit: log
    
    script:
    """
    trimmomatic PE -threads ${task.cpus} \
        ${reads[0]} ${reads[1]} \
        ${sample_id}_R1_paired.fq.gz ${sample_id}_R1_unpaired.fq.gz \
        ${sample_id}_R2_paired.fq.gz ${sample_id}_R2_unpaired.fq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
        2> ${sample_id}.trim.log
    """
}

process kraken2 {
    tag "${sample_id}"
    publishDir "${params.outdir}/kraken2", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path("${sample_id}.kraken2"), emit: classification
    path("${sample_id}.report"), emit: report
    
    script:
    def r1 = reads[0]
    def r2 = reads[1]
    """
    kraken2 --db ${params.kraken2_db} \
        --paired $r1 $r2 \
        --threads ${task.cpus} \
        --output ${sample_id}.kraken2 \
        --report ${sample_id}.report
    """
}


process spades_assembly {
    tag "${sample_id}"
    publishDir "${params.outdir}/assemblies", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.fasta"), emit: assembly

    script:
    def r1 = reads[0]
    def r2 = reads[1]
   
   """
    spades.py --careful \\
        -1 $r1 -2 $r2 \\
        -o ${sample_id}_assembly
	
	// To change fasta headers for downstream analysis
	
	if [[ -s ${sample_id}_assembly/scaffolds.fasta ]]; then
        # Rename headers to short format only if output is valid
        awk '/^>/ {print ">contig_" ++i} !/^>/ {print}' \\
            ${sample_id}_assembly/scaffolds.fasta > ${sample_id}.fasta
    else
        echo "ERROR: SPAdes failed to generate scaffolds.fasta"
        exit 1
    fi
    """
}

process quast {
    tag "${sample_id}"
    publishDir "${params.outdir}/quast", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    path("quast_results/report.html"), emit: report

    script:
    """
    quast.py -o quast_results $assembly
    """
}

process prokka {
    tag "${sample_id}"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("prokka_output/*"), emit: annotation

    conda '/home/anaconda3/envs/prokka'

    script:
    """
    prokka --outdir prokka_output \\
           --prefix ${sample_id} \\
           --cpus ${task.cpus} \\
           --force \\
           ${sample_id}.fasta
    """
}

process abricate {
    tag "${sample_id}"
    publishDir "${params.outdir}/abricate", mode: 'copy'
    
    input:
    tuple val(sample_id), path(assembly)
    
    output:
    path("*.tab"), emit: results
    
	// activate conda environment
	conda '/home/anaconda3/envs/abricate'
	
    script:
    """
    abricate --db ncbi $assembly > resistance_genes.tab
    abricate --db plasmidfinder $assembly > plasmids.tab
    """
}

// Workflow definition
workflow {
    // Quality control
    fastqc(reads)
    all_reports = fastqc.out.html.mix(fastqc.out.zip).collect()
    multiqc(all_reports)
    
    // Read trimming
    trimmomatic(reads)
    
    // Assembly
    spades_assembly(trimmomatic.out.trimmed_reads)
    
    // Run QUAST and Prokka in parallel
    //quast(spades_assembly.out.assembly)
	prokka(spades_assembly.out.assembly)
    //prokka(spades_assembly.out.assembly.map{ assembly -> [assembly.baseName.replaceFirst("_assembly/scaffolds",""), assembly] })
	abricate(spades_assembly.out.assembly)
}

// Customize the pipeline completion message
workflow.onComplete {
    println """
    Pipeline completed successfully!
    
    Results directory: ${params.outdir}
    Assemblies: ${params.outdir}/assemblies
    Annotations: ${params.outdir}/annotations
    
    To view quality reports:
    open ${params.outdir}/multiqc/multiqc_report.html
    open ${params.outdir}/quast/quast_results/report.html
    """
}
