# C.botulinum-C-gene-for-Neurtoxin

## Project 1 - Finding C.botulinum in Ancient DNA


# Quality Filtering and Trimming, Mapping Reads, and Generating Consensus Sequence
This repository contains a Bash script to perform quality filtering and trimming of sequencing data using FASTP, map reads to a reference genome using BWA, convert and manipulate files with SAMtools, and generate a consensus sequence using iVar. Additionally, it calculates the coverage of mapped reads and assesses the percent identity/similarity of the consensus sequence to the reference genome using BLAST.

## Prerequisites
Make sure you have the following tools installed:

fastp
bwa
samtools
ivar
blastn
You can install these tools via conda or your preferred package manager.

#Output Files
mapped_neurotoxin_data.sam: SAM file containing all mapped reads.
mapped_neurotoxin_data.bam: BAM file converted from the SAM file.
mapped_neurotoxin_data.sorted.bam: Sorted BAM file.
consensus_neurotoxin_sequence.fa: Consensus sequence generated by iVar.
results.txt: BLAST results showing percent identity/similarity.
output.text: File containing the coverage and BLAST results.

License

#Acknowledgments
This script and repository were created to facilitate the analysis of sequencing data, including quality control, mapping, and consensus sequence generation. Special thanks to the developers of FASTP, BWA, SAMtools, iVar, and BLAST for their powerful tools.

Feel free to modify and expand this script to suit your specific needs and datasets. If you encounter any issues or have suggestions for improvements, please open an issue or submit a pull request.
