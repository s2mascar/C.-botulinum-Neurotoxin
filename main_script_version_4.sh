#!/bin/bash
#SBATCH --job-name=neurotoxin_run
#SBATCH --account=def-acdoxey # adjust this to match the accounting group you are using to submit jobs
#SBATCH --time=5:00:00      # adjust this to match the walltime of your job
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --gres=gpu:0           # Request for GPUs
#SBATCH --cpus-per-task=8      # Request for CPUs
#SBATCH --mem=500G              # adjust this according to the memory requirement per node you need
#SBATCH --mail-user=s2mascar@uwaterloo.ca # adjust this to match your email address
#SBATCH --mail-type=ALL


echo "Clostridium Botulinum Neurotoxin Finder for Ancient DNA"

#To start the workflow, the respective SRA dataset needs to be downloaded and either seperated/not and also have the rest of the directories pre-set up

## HOW THIS WORKFLOW WORKS;
# 1. This script will be run from the scripts directory, it will need to access;
#		a)  Toxin reference seqeunce, SRA files
# 2. To start, we need to be in the directory of the SRA file, as it will have the seperate toxin gene directories for the data to be stored.
# 3. Once in the SRA directory, go to the "origin_files" folder and count how many files are present in the folder. (1 - Single-end, 2 - Paired-end)
# 4. In this "origin_files" folder, this is where the FASTP filtering will take place. Output should be SRA#_1_filtered.fastq etc... -- Maybe this should be its own step not in this script
# 5. Once FASTP is run, each toxin should have its own script, this script should be able to:
#		a) Find the toxin reference sequence, pointing to the same directory always (Folder: all_toxin_ref_seq)
#		b) Find the SRA files (Single or Paired), pointing to the "origin_files" (Folder: all_SRA_datasets) -- this is also the place where all the other mapped/BAM/SAM files will be housed.
#		c) The final output file (CSV) should contain the SRA name, toxin gene and also saved to (Folder: csv_files)
#		


#####################################################################################################################################

# In the SRA_datasets folder this is where all of the original downloaded SRA dataset should be stored, only place where they are stored and called from.
# In each SRA dataset folder, you should see how many files are there. 1 File - Single-End, 2 Files - Paired-End

# Define the directory
DIR="/home/smascar/scratch/neurotoxin_project/all_SRA_datasets/ERR5647172/origin_files"
SRA_name="ERR5647172"
toxin_ref_seq="/home/smascar/scratch/neurotoxin_project/all_toxin_ref_seq/C_reference_genome.fasta"
toxin_type="C"


# Count the number of regular files in the directory
#file_count=$(find "$DIR" -maxdepth 1 -type f | wc -l)

file_name_single=$(find "$SRA_name"_filtered.fastq)
file_name_paired=$(find "$SRA_name"_2_filtered.fastq)

# Check the number of files and take action accordingly
if [ "$file_name_single" -eq "$SRA_name"_filtered.fastq ]; then
	echo "This dataset is Single-End Reads"

	#################################################################################################

	#echo "Quality Filtering & Trimming via FASTP for Single-End Data"

	#fastp -i "$SRA_name".fastq -o "$SRA_name"_filtered.fastq


	echo "Mapping Single-End Reads to a Reference Genome"

	#Index the Reference Genome--------------------------------------------------------

	bwa index "$toxin_ref_seq"

	#Map Single-End Reads to the Reference Genome--------------------------------------------------

	bwa aln "$toxin_ref_seq" "$DIR"/"$SRA_name".fastq > "$SRA_name".sai

	#Making a .SAM file to contain all info about the Single-End read maps onto the ref genome--------------------

	bwa samse "$toxin_ref_seq" "$SRA_name".sai /data/"$SRA_name".fastq.gz > "$SRA_name".sam


elif [ "$file_name_paired" -eq "$SRA_name"_2_filtered.fastq ]; then
    echo "This dataset is Paired-End Reads"

    ##############################################################################################

    echo "Quality Filtering & Trimming via FASTP for Paired-End Data"

	fastp -i "$SRA_name"_1.fastq -I "$SRA_name"_2.fastq -o "$SRA_name"_1_filtered.fastq -O "$SRA_name"_2_filtered.fastq

	echo "Mapping Reads to a Reference Genome"

	#Index the Reference Genome--------------------------------------------------------

	bwa index C_reference_genome.fasta

	#Map Paired-End Reads to the Reference Genome--------------------------------------------------


	bwa aln C_reference_genome.fasta /data/Shyan/Project_2_Clostridium_botulinum_Neurotoxin/Full_Neurotoxin_Dataset_Analysis/SRA_datasets/"$SRA_name"_1_filtered.fastq > "$SRA_name"_1_filtered.sai

	bwa aln C_reference_genome.fasta /data/Shyan/Project_2_Clostridium_botulinum_Neurotoxin/Full_Neurotoxin_Dataset_Analysis/SRA_datasets/"$SRA_name"_2_filtered.fastq > "$SRA_name"_2_filtered.sai


	#Making a .SAM file to contain all info about the Paired-End read maps onto the ref genome--------------------

	
	bwa sampe C_reference_genome.fasta "$SRA_name"_1_filtered.sai "$SRA_name"_2_filtered.sai /data/Shyan/Project_2_Clostridium_botulinum_Neurotoxin/Full_Neurotoxin_Dataset_Analysis/SRA_datasets/"$SRA_name"_1_filtered.fastq /data/Shyan/Project_2_Clostridium_botulinum_Neurotoxin/Full_Neurotoxin_Dataset_Analysis/SRA_datasets/"$SRA_name"_2_filtered.fastq > mapped_neurotoxin_data.sam

else
    echo "There are no files to run. Error."
fi



######################################### LEAVE THIS PART IT DOESN'T REQUIRE WHETHER SINGLE OR PAIRED END


#Index the Reference genome to work with Samtools

samtools faidx C_reference_genome.fasta


#Convert the .SAM file to a .BAM file

samtools view -S -b mapped_neurotoxin_data.sam > mapped_neurotoxin_data.bam


#Sort the .BAM file and index it
samtools sort mapped_neurotoxin_data.bam > mapped_neurotoxin_data.sorted.bam
samtools index mapped_neurotoxin_data.sorted.bam



echo "Coverage of Mapped Reads to Reference Genome" >> "$SRA_name"_"$toxin_type"_output.txt

samtools coverage mapped_neurotoxin_data.sorted.bam >> "$SRA_name"_"$toxin_type"_output.txt


echo "Make a Consensus Sequence"

samtools mpileup -aa -A -Q 0 -d 0 mapped_neurotoxin_data.sorted.bam | ivar consensus -p consensus_neurotoxin_sequence -m 2 -n N -t 0.5


echo "Percent Identity/Similarity from Consensus Sequence to Reference Genome" >> "$SRA_name"_"$toxin_type"_output.txt

#Check Python & Biopython is installed as well.

python3 --version

pip3 install biopython

python3 calculateIdentity.py aligned_sequences.fasta >> "$SRA_name"_"$toxin_type"_output.txt


echo "Extract the necessary data into a CSV and move the CSV file to the all_CSV_files directory"

pip3 install pandas

python3 get_coverage_value.py "$SRA_name"_"$toxin_type"_output.txt

DESTINATION_DIR="/home/smascar/scratch/neurotoxin_project/csv_files"

# Copy the output file to the destination directory

cp "$SRA_name"_"$toxin_type"_output.txt "$DESTINATION_DIR" 







