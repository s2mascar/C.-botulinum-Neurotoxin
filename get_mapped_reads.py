#!/usr/bin/env python3

import os
import subprocess

def execute_bash_command_in_directories(root_dir, target_file_name, command_template):
    """
    Recursively execute a bash command in directories containing a specific file.
    
    Parameters:
    - root_dir: str: The root directory to start searching from.
    - target_file_name: str: The name of the file to look for.
    - command_template: str: The bash command template to execute.
    """
    for dirpath, dirnames, filenames in os.walk(root_dir):
        if target_file_name in filenames:
            file_path = os.path.join(dirpath, target_file_name)
            command = command_template.format(file_path=file_path)
            try:
                # Change directory to dirpath before executing the command
                subprocess.run(command, shell=True, check=True, cwd=dirpath)
                print(f"Executed in {dirpath}: {command}")
            except subprocess.CalledProcessError as e:
                print(f"Error executing command in {dirpath}: {e}")

if __name__ == "__main__":
    # Define the root directory, target file, and command template
    root_directory = "/home/smascar/scratch/neurotoxin_project/all_SRA_datasets"
    target_file = "mapped_neurotoxin_data.sorted.bam"
    bash_command_template = "samtools view -b -F 4 {file_path} > mapped.bam"
    
    # Call the function to execute the command
    execute_bash_command_in_directories(root_directory, target_file, bash_command_template)





#Extract the mapped reads from the sorted BAM file
  samtools view -b -F 4 mapped_neurotoxin_data.sorted.bam > only__mappedreads_neurotoxin_sorted.bam

root_directory = "/home/smascar/scratch/neurotoxin_project/all_SRA_datasets/ERR1879296/Type_A"
    