#!/usr/bin/bash

echo "Clostridium Botulinum Neurotoxin Finder for Ancient DNA"

#To start the workflow, make the following directories for the data to be stored into.

mkdir ERR5647172
cd ERR5647172
mkdir Type_A Type_B Type_C Type_D Type_E Type_F Other_types
cd Type_A
mkdir A1 A2 A3 A4 A5 A6 A7 A8

cd Type_B
mkdir B1 B2 B3 B4 B5 B6 B7 B8

cd Type_C
mkdir C CD

cd Type_D
mkdir D DC

cd Type_E
mkdir E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12

cd Type_F
mkdir F1 F2 F3 F4 F5 F6 F7 F8 F9

cd Other_types
mkdir En G HA_FA HDR7951433 PMP1 TeNT Wo X


# In the SRA_datasets folder this is where all of the original downloaded SRA dataset should be stored, only place where they are stored and called from.
# In each SRA dataset folder, you should see how many files are there. 1 File - Single-End, 2 Files - Paired-End

# Define the directory
DIR="/data/Shyan/Project_2_Clostridium_botulinum_Neurotoxin/Full_Neurotoxin_Dataset_Analysis/SRA_datasets/ERR5647172"

# Count the number of regular files in the directory
file_count=$(find "$DIR" -maxdepth 1 -type f | wc -l)

# Check the number of files and take action accordingly
if [ "$file_count" -eq 1 ]; then
	echo "This dataset is Single-End Reads"

	#################################################################################################

	echo "Quality Filtering & Trimming via FASTP for Single-End Data"

	fastp -i ERR5647172.fastq -o ERR5647172_filtered.fastq


	echo "Mapping Single-End Reads to a Reference Genome"

	#Index the Reference Genome--------------------------------------------------------

	bwa index F1_reference_genome.fasta

	#Map Single-End Reads to the Reference Genome--------------------------------------------------

	bwa aln ecoli-rel606.fa /data/SRR098038.fastq.gz > SRR098038.sai

	#Making a .SAM file to contain all info about the Single-End read maps onto the ref genome--------------------

	bwa samse ecoli-rel606.fa SRR098038.sai /data/SRR098038.fastq.gz > SRR098038.sam




elif [ "$file_count" -eq 2 ]; then
    echo "This dataset is Paired-End Reads"

    ##############################################################################################

    echo "Quality Filtering & Trimming via FASTP for Paired-End Data"

	fastp -i ERR5647172_1.fastq -I ERR5647172_2.fastq -o ERR5647172_1_filtered.fastq -O ERR5647172_2_filtered.fastq

	echo "Mapping Reads to a Reference Genome"

	#Index the Reference Genome--------------------------------------------------------

	bwa index F1_reference_genome.fasta

	#Map Paired-End Reads to the Reference Genome--------------------------------------------------


	bwa aln F1_reference_genome.fasta /data/Shyan/Project_2_Clostridium_botulinum_Neurotoxin/Full_Neurotoxin_Dataset_Analysis/SRA_datasets/ERR5647172_1_filtered.fastq > ERR5647172_1_filtered.sai

	bwa aln F1_reference_genome.fasta /data/Shyan/Project_2_Clostridium_botulinum_Neurotoxin/Full_Neurotoxin_Dataset_Analysis/SRA_datasets/ERR5647172_2_filtered.fastq > ERR5647172_2_filtered.sai


	#Making a .SAM file to contain all info about the Paired-End read maps onto the ref genome--------------------

	
	bwa sampe F1_reference_genome.fasta ERR5647172_1_filtered.sai ERR5647172_2_filtered.sai /data/Shyan/Project_2_Clostridium_botulinum_Neurotoxin/Full_Neurotoxin_Dataset_Analysis/SRA_datasets/ERR5647172_1_filtered.fastq /data/Shyan/Project_2_Clostridium_botulinum_Neurotoxin/Full_Neurotoxin_Dataset_Analysis/SRA_datasets/ERR5647172_2_filtered.fastq > mapped_neurotoxin_data.sam

else
    echo "There are no files to run. Error."
fi

#__________________________________________________________________________________________________________________________________________

echo "Quality Filtering & Trimming via FASTP"

#Single-End Data

#fastp -i in.fq -o out.fq

#Paired-End Data

#fastp -i ERR5647172_1.fastq -I ERR5647172_2.fastq -o ERR5647172_1_filtered.fastq -O ERR5647172_2_filtered.fastq



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




######################################### LEAVE THIS PART IT DOESN'T REQUIRE WHETHER SINGLE OR PAIRED END


#Index the Reference genome to work with Samtools

samtools faidx F1_reference_genome.fasta


#Convert the .SAM file to a .BAM file

samtools view -S -b mapped_neurotoxin_data.sam > mapped_neurotoxin_data.bam


#Sort the .BAM file and index it
samtools sort mapped_neurotoxin_data.bam > mapped_neurotoxin_data.sorted.bam
samtools index mapped_neurotoxin_data.sorted.bam



echo "Coverage of Mapped Reads to Reference Genome" >> output.txt

samtools coverage mapped_neurotoxin_data.sorted.bam >> output.txt


echo "Make a Consensus Sequence"

samtools mpileup -aa -A -Q 0 -d 0 mapped_neurotoxin_data.sorted.bam | ivar consensus -p consensus_neurotoxin_sequence -m 2 -n N -t 0.5


echo "Percent Identity/Similarity from Consensus Sequence to Reference Genome" >> output.text

#Check Python & Biopython is installed as well.

python3 --version

pip3 install biopython

python3 calculateIdentity.py aligned_sequences.fasta >> output.txt


echo "Extract the necessary data into a CSV and move the CSV file to the all_CSV_files directory"

pip3 install pandas

python3 get_coverage_value.py output.text

DESTINATION_DIR="/data/Shyan/Project_2_Clostridium_botulinum_Neurotoxin/Full_Neurotoxin_Dataset_Analysis"

# Copy the output file to the destination directory

cp output.txt "$DESTINATION_DIR" 







