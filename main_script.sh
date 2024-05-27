#!/usr/bin/bash

echo "Quality Filtering & Trimming via FASTP"

#Single-End Data

#fastp -i in.fq -o out.fq

#Paired-End Data

#fastp -i ERR5647172_1_filtered.fastq -I ERR5647172_2_filtered.fastq -o ERR5647172_1_filtered.fastq -O ERR5647172_2_filtered.fastq



echo "Mapping Reads to a Reference Genome"

#Index the Reference Genome--------------------------------------------------------

bwa index F1_reference_genome.fasta

#Map Reads to the Reference Genome--------------------------------------------------

#Single-End Data

#bwa aln ecoli-rel606.fa /data/SRR098038.fastq.gz > SRR098038.sai


#Paired-End Data

bwa aln F1_reference_genome.fasta /data/Shyan/Project_2_Clostridium_botulinum_Neurotoxin/Full_Neurotoxin_Dataset_Analysis/SRA_datasets/ERR5647172_1_filtered.fastq > ERR5647172_1_filtered.sai

bwa aln F1_reference_genome.fasta /data/Shyan/Project_2_Clostridium_botulinum_Neurotoxin/Full_Neurotoxin_Dataset_Analysis/SRA_datasets/ERR5647172_2_filtered.fastq > ERR5647172_2_filtered.sai

#/data/Shyan/Project_2_Clostridium_botulinum_Neurotoxin/Full_Neurotoxin_Dataset_Analysis/SRA_datasets

#Making a .SAM file to contain all info about read maps onto the ref genome--------------------

#Single-End Data
#bwa samse ecoli-rel606.fa SRR098038.sai /data/SRR098038.fastq.gz > SRR098038.sam

#Paired-End Data
bwa sampe F1_reference_genome.fasta ERR5647172_1_filtered.sai ERR5647172_2_filtered.sai /data/Shyan/Project_2_Clostridium_botulinum_Neurotoxin/Full_Neurotoxin_Dataset_Analysis/SRA_datasets/ERR5647172_1_filtered.fastq /data/Shyan/Project_2_Clostridium_botulinum_Neurotoxin/Full_Neurotoxin_Dataset_Analysis/SRA_datasets/ERR5647172_2_filtered.fastq > mapped_neurotoxin_data.sam


#Index the Reference genome to work with Samtools

samtools faidx F1_reference_genome.fasta


#Convert the .SAM file to a .BAM file

samtools view -S -b mapped_neurotoxin_data.sam > mapped_neurotoxin_data.bam


#Sort the .BAM file and index it
samtools sort mapped_neurotoxin_data.bam > mapped_neurotoxin_data.sorted.bam
samtools index mapped_neurotoxin_data.sorted.bam



echo "Coverage of Mapped Reads to Reference Genome" >> output.text

samtools coverage mapped_neurotoxin_data.sorted.bam >> output.text


echo "Make a Consensus Sequence"

samtools mpileup -aa -A -Q 0 -d 0 mapped_neurotoxin_data.sorted.bam | ivar consensus -p consensus_neurotoxin_sequence -m 2 -n N -t 0.5


echo "Percent Identity/Similarity from Consensus Sequence to Reference Genome" >> output.text

blastn -query F1_reference_genome.fasta -subject consensus_neurotoxin_sequence.fa -outfmt 6 -out results.txt









