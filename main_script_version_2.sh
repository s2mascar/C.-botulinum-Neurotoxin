#!/usr/bin/bash

echo "Clostridium Botulinum Neurotoxin Finder for Ancient DNA"

#To start the workflow, download the respective SRA dataset

#use a loop for each of the SRA filenames


list_of_sra_datasets=(ERR5647172 ERR4375049 ERR4374027 SRR23016251 SRR23016252 SRR1748806 SRR9276195 SRR9276194 SRR11615772 ERR1879296 ERR3678614 ERR3678629 ERR3678626)


for dataset in list_of_sra_datasets; do
  # Commands to be executed for each item
  echo "Processing item: $item"
done




#####################################################################################################################################

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

DESTINATION_DIR="/data/Shyan/Project_2_Clostridium_botulinum_Neurotoxin/Full_Neurotoxin_Dataset_Analysis/all_CSV_files"

# Copy the output file to the destination directory

cp output.txt "$DESTINATION_DIR" 







