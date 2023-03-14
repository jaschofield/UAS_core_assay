# UAS_core_assay
Code associated with Schofield and Hahn 2023 publication.

Information about data:
* Transcriptional reporter plasmid libraries were generated with thousands of Upstream Activating Sequences (UASs) cloned upstream of 4 representative core promoters from TFIID and CR classes of yeast genes. UAS-core combinations drive expression of a reporter gene containing a random sequence barcode. Nascent transcription is measured through targeted nascent RNA-seq of reporter barcodes, and assays were performed with and without rapid depletion of transcriptional coactivator complexes (plus/minus auxin in Med15, Spt7, or Taf13 degron strains). To assign barcodes to their UAS-core combinations of origin and normalize to reporter DNA input, DNA-seq was performed on plasmid libraries purified from transformed yeast.

Information about the analysis:
* This code uses VSEARCH (torognes) clustering to identify and count consensus barcode and UAS-core sequences, allowing sequences with potential sequencing mutations to be grouped into appropriate consensus clusters. The clustering is performed in sequential steps to reduce computational time, and the strictness of clustering is higher for barcode sequences to avoid combining expression from highly similar barcodes. A fuzzy join is used to identify the core promoter and UAS associated with each DNA consensus cluster, and an exact join is used to match RNA barcodes with their DNA consensus clusters. Ambiguous barcodes (those assigned to more than one UAS-core combination) are removed from analysis, and expression levels for each UAS-core combination are determined by normalizing RNA barcode counts to library size and DNA input.

Instructions for running code:
1. Ensure VSEARCH is installed.
2. Ensure the following R packages and their dependencies are installed: tidyverse, purrr, Biostrings, fuzzyjoin.
3. Download FastQ files from GEO repository (GSE217230).
4. Run RNA_barcode_format.sh for each RNA-seq sample.
5. Run DNA_fa_join_ex.sh to execute DNA_fa_join.R script for each DNA-seq sample.
6. Run Fasta_split_RNA.sh for outputs from RNA_barcode_format, and Fasta_split_DNA.sh for outputs from DNA_fa_join.
7. Run Cluster_1_RNA.sh for each split from Fasta_split_RNA, and Cluster_1_DNA.sh for each split from Fasta_split_DNA.
* Split names are variable for Cluster 1 scripts. An example of running script: sbatch -c N Cluster_1_DNA.sh Sample_split_aa.
* Record name of each slurm.out file produced for each sample split. Lists of slurm.out file names for each sample will be input into next step.
8. Run Cluster_2_RNA.sh and Cluster_2_DNA.sh using lists of slurm.out for each sample.
9. Run Master_R_process_ex.sh to execute Master_R_process.R to perform final processing of data.
