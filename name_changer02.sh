#!/bin/bash

# Define the directory to search
DIRECTORY_1="/data/Shyan/trial"

# Define the text to search for and replace
OLD_TEXT_1=""
NEW_TEXT_2="B2_reference_genome"

# Find all files and apply the sed command
find "$DIRECTORY_1" -type f -name '*.sh' -exec sed -i "s/$OLD_TEXT_1/$NEW_TEXT_1/g" {} +

echo "Replacement 1 complete."

# Define the directory to search
DIRECTORY_2="/data/Shyan/trial"

# Define the text to search for and replace
OLD_TEXT_2="SRAname_b2_output"
NEW_TEXT_2="SRAname_b2_SRAname_b2_output"

# Find all files and apply the sed command
find "$DIRECTORY_2" -type f -name '*.sh' -exec sed -i "s/$OLD_TEXT_2/$NEW_TEXT_2/g" {} +

echo "Replacement 2 complete."

list_of_sra_datasets=(ERR5647172 ERR4375049 ERR4374027 SRR23016251 SRR23016252 SRR1748806 SRR9276195 SRR9276194 SRR11615772 ERR1879296 ERR3678614 ERR3678629 ERR3678626)

