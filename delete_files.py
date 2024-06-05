import os

def delete_specific_file(root_dir, file_name_to_delete):
    """
    Recursively delete specific files in a directory and its subdirectories.
    
    Parameters:
    - root_dir: str: The root directory to start searching from.
    - file_name_to_delete: str: The name of the file to delete.
    """
    for dirpath, dirnames, filenames in os.walk(root_dir):
        for file in filenames:
            if file == file_name_to_delete:
                file_path = os.path.join(dirpath, file)
                try:
                    os.remove(file_path)
                    print(f"Deleted: {file_path}")
                except Exception as e:
                    print(f"Error deleting {file_path}: {e}")

if __name__ == "__main__":
    # Define the root directory and the file name to delete
    root_directory = "/data/Shyan/trial"
    file_name = "ex1"
    
    # Call the function to delete the specific file
    delete_specific_file(root_directory, file_name)
