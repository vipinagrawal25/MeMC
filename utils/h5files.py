import os
import shutil
from pathlib import Path
import glob
import sys

def get_all_folders(directory):
    folder_names = []
    # Traverse the directory recursively
    for root, dirs, files in os.walk(directory):
        # Check if 'mc_log' file is in the current folder and 'cs0o10' is in the path
        if 'mc_log' in files:
            folder_names.append(root)
    return folder_names

def lastfile(dirname, prefix='', zfill=5, suffix=".h5"):
    dirname = os.path.abspath(dirname)
    if not dirname.endswith(os.sep):
        dirname += os.sep

    # Find all files matching the pattern
    file_pattern = f"{dirname}{prefix}*{suffix}"
    files = glob.glob(file_pattern)

    # If no files found, return None
    if not files:
        return None
    try:
        files_sorted = sorted(
            files,
            key=lambda f: int(os.path.basename(f).replace(prefix, '').replace(suffix, ''))
        )
        return files_sorted[-1]
    except ValueError:
        return None

source_root = sys.argv[1]  
target_folder = source_root+"/latest_files"  

# Create the target folder if it doesn't exist
os.makedirs(target_folder, exist_ok=True)
alldirs = get_all_folders(source_root)
# print(alldirs)

for dd in alldirs:
    latest_h5_path=lastfile(dd,"snap_")
    # Create a target name based only on the relative path to the folder
    relative_path = os.path.relpath(dd, source_root)
    sanitized_path = relative_path.replace(os.sep, '_').replace('.', 'o')
    target_h5_name = f"{sanitized_path}.h5"
    target_h5_path = os.path.join(target_folder, target_h5_name)
    print(latest_h5_path, target_h5_path)
    shutil.copy(latest_h5_path, target_h5_path)