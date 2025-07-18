Nextflow pipeline for bacterial whole genome sequence analysis </h2>

<h4>Software included in the pipeline.</h4>

<h4>Tools used</h4>
Quality Control<br>
a) <h4></h4>FastQC</h4> <br> is widely used and is robust, efficient, and versatile quality control software for a varied range of raw genetic data. It outputs a quality report which can be viewed to give an indication on how good the respective reads are. For more information, please visit the FastQC website.

b) <h4>MultiQC</h4> <br> is used create a single report with interactive plots for multiple bioinformatics analyses across many samples. It reports can describe multiple analysis steps and large numbers of samples within a single plot, and multiple analysis tools making it ideal for routine fast quality control. For more information, please visit the MultiQC website.

c) <h4>Trimmomatic</h4> is a commonly used tool to remove adapter sequences, low-quality bases, and low-quality reads from sequencing reads. For more information and usage instructions please visit the Trimmomatic website.

<h4>SPAdes</h4>
It was used to perform de novo genome assembly after quality control and trimming the adapters from the raw reads. SPAdes is a de Bruijn graph based assembler, designed and intended for small genomes and can take sequencing data from Illumina and IonTorrent platforms. For more information, please visit the Spades website.

<h4>Prokka</h4>
It is a software tool to annotate bacterial, archaeal and viral genomes quickly and produce standards-compliant output files. It identifies features of interest in a set of genomic DNA sequences, and labelling them with useful information. For more information, please visit the Prokka website.

<h4>QUAST</h4>
It is a quality assessment tool for evaluating and comparing genome assemblies by computing various metrics and works both with and without reference genomes. It produces many reports, summary tables and plots to help scientists in their research and in their publications. For more information, please visit the Quast website.

<h4>Abricate</h4>
It is a tool for mass screening of contigs for antimicrobial resistance or virulence genes. It comes bundled with multiple databases: NCBI, CARD, ARG-ANNOT, Resfinder, MEGARES, EcOH, PlasmidFinder, Ecoli_VF and VFDB. For more information, please visit https://github.com/tseemann/abricate.


