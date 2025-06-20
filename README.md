<h2> Nextflow pipeline for bacterial whole genome sequence analysis </h2>

<h3>Quality control fo the raw sequencing data generated from the sequencing instrument to remove low-quality reads and adapter sequences.</h3>

1). FastQC is used to assess the quality distribution and other quality metrics of sequencing reads. It helps detect potential issues during the sequencing process, such as decreasing sequencing quality and adapter contamination. For more information, please visit the FastQC website.
2). MultiQC is a tool to create a single report with interactive plots for multiple bioinformatics analyses across many samples. MultiQC reports can describe multiple analysis steps and large numbers of samples within a single plot, and multiple analysis tools making it ideal for routine fast quality control.For more information, please visit the MultiQC website.
3). Trimmomatic: Used to remove adapter sequences, low-quality bases, and low-quality reads from sequencing reads. You can find more information and usage instructions on the Trimmomatic GitHub.
4). Kraken 2 is the newest version of Kraken, a taxonomic classification system using exact k-mer matches to achieve high accuracy and fast classification speeds. This classifier matches each k-mer within a query sequence to the lowest common ancestor (LCA) of all genomes containing the given k-mer. For more information, please visit the Kraken 2 website. 
