







list_of_sra_datasets=(ERR4375049 ERR4374027 SRR23016251 SRR23016252 SRR1748806 SRR9276195 SRR9276194 SRR11615772 ERR1879296 ERR3678614 ERR3678629 ERR3678626)


list of 
list_of_toxin_subtype=(A1 A2 A3 A4 A5 A6 A7 A8 B1 B2 B3 B4 B5 B6 B7 B8)
DIR="/home/smascar/scratch/neurotoxin_project/all_SRA_datasets"

SRA_name="ERR5647172"

toxin_C_DIR="Type_C/CD"
toxin_type="CD"

ref_seq="CD_reference_genome.fasta"


for dataset in "${list_of_sra_datasets[@]}"; do

	cd "$DIR"/"$dataset"/"$toxin_C_DIR"
	
	samtools view -b -F 4 file.bam > mapped.bam



done
