<h1> Nextflow pipeline for bacterial whole genome sequence analysis </h1>

<h2>Quality control fo the raw sequencing data generated from the sequencing instrument to remove low-quality reads and adapter sequences.</h2>

    Trimmomatic: Used to remove adapter sequences, low-quality bases, and low-quality reads from sequencing reads. You can find more information and usage instructions on the Trimmomatic GitHub.

    FastQC: Used to assess the quality distribution and other quality metrics of sequencing reads. It helps detect potential issues during the sequencing process, such as decreasing sequencing quality and adapter contamination. For more information, please visit the FastQC website.

    BWA-MEM: Used to align sequencing reads to a reference genome. BWA-MEM employs a fast seed-and-extend algorithm that efficiently handles longer reads. Detailed documentation and usage instructions can be found on the BWA-MEM GitHub.

    Bowtie2: Used to align reads to a reference genome. It combines a fast seed-and-extend algorithm with a maximum-likelihood alignment algorithm, suitable for aligning short sequences. More information can be found on the Bowtie2 website.

    SAMtools: A toolkit for processing alignment results, which can perform operations such as conversion, sorting, indexing, and filtering of SAM/BAM format files. Additionally, it provides a set of command-line tools for statistics on alignment information, SNP/Indel detection, and more. More documentation and download links can be found on the SAMtools website.

    GATK (Genome Analysis Toolkit): A toolkit for analyzing and detecting variants in genomic data. It offers various algorithms and tools, including genome re-alignment, SNP/Indel detection, and variant filtering. For detailed information, please visit the GATK website.

    Picard Tools: A toolkit for processing SAM/BAM files, which can mark PCR duplicates, calculate coverage, create statistical reports, and more. Additional information and usage instructions can be found on the Picard Tools website.

    QualiMap: A tool used to assess the quality of sequencing data. It detects quality metrics such as coverage, GC content, error rates of sequencing reads, and generates corresponding statistical charts and reports. For more information, please visit the QualiMap website.
