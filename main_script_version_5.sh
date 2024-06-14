#!/bin/bash
#SBATCH --job-name=neurotoxin_run
#SBATCH --account=def-acdoxey # adjust this to match the accounting group you are using to submit jobs
#SBATCH --time=5:00:00      # adjust this to match the walltime of your job
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4      # Request for CPUs
#SBATCH --mem=500G              # adjust this according to the memory requirement per node you need
#SBATCH --mail-user=s2mascar@uwaterloo.ca # adjust this to match your email address
#SBATCH --mail-type=ALL

echo "Toxin En Script -- ALL DATASETS"
echo "Load the Modules"

module load StdEnv/2023 bwa samtools ivar

DIR="/home/smascar/scratch/neurotoxin_project/all_SRA_datasets"

list_of_sra_datasets=(ERR4375049 ERR4374027 SRR23016251 SRR23016252 SRR1748806 SRR9276195 SRR9276194 SRR11615772 ERR1879296 ERR3678614 ERR3678629 ERR3678626 ERR5647172)

toxin_C_DIR="Other_types/En"
toxin_type="En"

ref_seq="En_reference_genome.fasta"


for dataset in "${list_of_sra_datasets[@]}"; do
  # Commands to be executed for each item
  echo "Processing item: $dataset"
  
  #echo "Quality Filtering & Trimming via FASTP for Paired-End Data"

  #cd "$DIR"/"$dataset"
  
  #fastp -i "$dataset"_1.fastq -I "$dataset"_2.fastq -o "$dataset"_1_filtered.fastq -O "$dataset"_2_filtered.fastq


  echo "Mapping Reads to a Reference Genome"

  cd "$DIR"/"$dataset"/"$toxin_C_DIR"
  #Index the Reference Genome--------------------------------------------------------

  bwa index /home/smascar/scratch/neurotoxin_project/all_toxin_ref_seq/"$ref_seq"

  #Map Paired-End Reads to the Reference Genome--------------------------------------------------

  bwa aln /home/smascar/scratch/neurotoxin_project/all_toxin_ref_seq/"$ref_seq" "$DIR"/"$dataset"/"$dataset"_1_filtered.fastq > "$dataset"_1_filtered.sai

  bwa aln /home/smascar/scratch/neurotoxin_project/all_toxin_ref_seq/"$ref_seq" "$DIR"/"$dataset"/"$dataset"_2_filtered.fastq > "$dataset"_2_filtered.sai


  #Making a .SAM file to contain all info about the Paired-End read maps onto the ref genome--------------------


  bwa sampe /home/smascar/scratch/neurotoxin_project/all_toxin_ref_seq/"$ref_seq" "$dataset"_1_filtered.sai "$dataset"_2_filtered.sai "$DIR"/"$dataset"/"$dataset"_1_filtered.fastq "$DIR"/"$dataset"/"$dataset"_2_filtered.fastq > mapped_neurotoxin_data.sam


  #Index the Reference genome to work with Samtools

  samtools faidx /home/smascar/scratch/neurotoxin_project/all_toxin_ref_seq/"$ref_seq"


  #Convert the .SAM file to a .BAM file


  samtools view -S -b mapped_neurotoxin_data.sam > mapped_neurotoxin_data.bam


  #Sort the .BAM file and index it


  samtools sort mapped_neurotoxin_data.bam > mapped_neurotoxin_data.sorted.bam
  samtools index mapped_neurotoxin_data.sorted.bam



  echo "Coverage of Mapped Reads to Reference Genome" >> "$dataset"_"$toxin_type"_output.txt

  samtools coverage mapped_neurotoxin_data.sorted.bam >> "$dataset"_"$toxin_type"_output.txt


  echo "Make a Consensus Sequence - with only 1 Read matching"

  samtools mpileup -aa -A -Q 0 -d 0 mapped_neurotoxin_data.sorted.bam | ivar consensus -p consensus_neurotoxin_sequence -m 1 -n N -t 0.5

  #Extract the mapped reads from the sorted BAM file
  
  samtools view -b -F 4 mapped_neurotoxin_data.sorted.bam > only__mappedreads_neurotoxin_sorted.bam

  #Index the BAM file to use Samtools tview

  samtools index only__mappedreads_neurotoxin_sorted.bam
  ################################################################################################################

  echo "Percent Identity/Similarity from Consensus Sequence to Reference Genome" >> "$dataset"_"$toxin_type"_output.txt

  # Add the Consensus Sequence and the Reference Sequence Together

  cat /home/smascar/scratch/neurotoxin_project/all_toxin_ref_seq/"$ref_seq" consensus_neurotoxin_sequence.fa > aligned_sequences.fasta

  #Check Python & Biopython is installed as well.

  python3 --version

  pip3 install biopython

  python3 /home/smascar/scratch/neurotoxin_project/scripts/calculateIdentity.py aligned_sequences.fasta >> "$dataset"_"$toxin_type"_output.txt

  echo "Extract the necessary data into a CSV and move the CSV file to the all_CSV_files directory"

  pip3 install pandas

  python3 /home/smascar/scratch/neurotoxin_project/scripts/get_coverage_value.py "$dataset"_"$toxin_type"_output.txt

  DESTINATION_DIR="/home/smascar/scratch/neurotoxin_project/csv_files"

  # Copy the output file to the destination directory

  cp summary_output_"$dataset"_"$toxin_type"_output.txt.csv "$DESTINATION_DIR" 


  echo "Delete Large Files"

  #mapped_neurotoxin_data.bam, mapped_neurotoxin_data.sam, mapped_neurotoxin_data.sorted.bam

  rm mapped_neurotoxin_data.bam

  rm mapped_neurotoxin_data.sam

  rm mapped_neurotoxin_data.sorted.bam


done



