import os
import shutil
import f90nml

def get_value_from_namelist(folder, keys):
    """Reads the specified keys from the namelist in a folder."""
    namelist_path = os.path.join(folder, "para_file.in")
    if not os.path.exists(namelist_path):
        return None
    
    para = f90nml.read(namelist_path)
    return {key: para[section][param] for key, (section, param) in keys.items()}

def sanitize_folder_name(value):
    """Converts a float value to a string with 'o' replacing the decimal point."""
    return f"{value:.2f}".replace('.', 'o')

# Keys to extract from the namelist
keys_to_extract = {
    "charge": ("electrostatpara", "charge2"),
    "conc": ("electrostatpara", "conc"),
    "compfrac": ("meshpara", "compfrac"),
}

# Root directories
source_root = "./"  # Adjust if needed
target_root = "./"

# Create target root directory if it doesn't exist
os.makedirs(target_root, exist_ok=True)

# Iterate through all subdirectories
for folder in sorted(os.listdir(source_root)):
    source_folder = os.path.join(source_root, folder)
    if not os.path.isdir(source_folder):
        continue

    # Extract parameter values
    values = get_value_from_namelist(source_folder, keys_to_extract)
    if values is None:
        print(f"Skipping folder {source_folder}, no 'para_file.in' found.")
        continue

    # Generate the target folder path
    charge_dir = f"ch{sanitize_folder_name(values['charge'])}"
    conc_dir = f"cs{sanitize_folder_name(values['conc']*1000)}"
    compfrac_dir = f"fr{sanitize_folder_name(values['compfrac'])}"
    target_folder = os.path.join(target_root, charge_dir, conc_dir, compfrac_dir)

    # Create target directories if they don't exist
    os.makedirs(target_folder, exist_ok=True)

    # Move contents of the source folder to the target folder
    for item in os.listdir(source_folder):
        source_item_path = os.path.join(source_folder, item)
        target_item_path = os.path.join(target_folder, item)
        shutil.move(source_item_path, target_item_path)
    print(f"Moved contents of {source_folder} to {target_folder}.")