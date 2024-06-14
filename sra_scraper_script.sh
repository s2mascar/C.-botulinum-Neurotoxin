#!/usr/bin/bash

echo "Clostridium Botulinum Neurotoxin Finder for Ancient DNA"

#To start the workflow, download the respective SRA dataset

#use a loop for each of the SRA filenames

list_of_sra_datasets=(ERR5647172 ERR4375049 ERR4374027 SRR23016251 SRR23016252 SRR1748806 SRR9276195 SRR9276194 SRR11615772 ERR1879296 ERR3678614 ERR3678629 ERR3678626)

for dataset in "${list_of_sra_datasets[@]}"; do
  # Commands to be executed for each item
  echo "Processing item: $dataset"

  # Download the SRA file using prefetch
  prefetch "$dataset"

  # Convert to fastq using fasterq-dump with split-files option
  fasterq-dump "$dataset" --split-files
done




