#!/bin/bash

# Define the output filename
#OUTPUT_FILENAME="output.text"  # Replace this with the actual filename if needed

# Run the Python script with the output filename as an argument
#python3 get_coverage_value.py "$OUTPUT_FILENAME"

# Change to the target directory
cd /data/Shyan/Project_2_Clostridium_botulinum_Neurotoxin/Full_Neurotoxin_Dataset_Analysis/SRA_datasets/ERR5647172

# Define the directory
DIR="/data/Shyan/Project_2_Clostridium_botulinum_Neurotoxin/Full_Neurotoxin_Dataset_Analysis/SRA_datasets/ERR5647172"

# Count the number of regular files in the directory
file_count=$(find "$DIR" -maxdepth 1 -type f | wc -l)

# Check the number of files and take action accordingly
if [ "$file_count" -eq 1 ]; then
    echo "There is 1 file in the directory."
elif [ "$file_count" -eq 2 ]; then
    echo "There are 2 files in the directory."
else
    echo "There are $file_count files in the directory."
fi


