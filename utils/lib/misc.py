import numpy as np
import h5py
from numpy import linalg as LA
import os
import glob
import re
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
home_dir = os.path.expanduser("~")
plt.style.use(os.path.join(home_dir, '.matplotlibrc'))

def lastfilename(dirname, prefix='snap_', zfill=5, suffix=".vtk"):
    """
    Returns the full path of the last file matching the pattern in a directory.
    
    Parameters:
        - dirname: Directory path.
        - prefix: Prefix of the file name (default: '').
        - zfill: Zero-padded width of the numerical identifier (default: 5).
        - suffix: File extension (default: '.h5').
    
    Returns:
        - str: Full path of the last file.
        - None: If no files are found.
    """
    # Normalize the directory path
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
        # Extract numeric parts of the file names and sort by them
        files_sorted = sorted(
            files,
            key=lambda f: int(os.path.basename(f).replace(prefix, '').replace(suffix, ''))
        )
        # Return the last file in the sorted list
        return files_sorted[-1]
    except ValueError:
        # If parsing fails, return None
        return None

def extract_parameters(filename):
    """
    Extract charge, concentration, and fraction from the filename.
    Example format: ch0o03_cs0o01_fr0o50.png

    Parameters:
        filename (str): Filename to parse.

    Returns:
        tuple: (charge, concentration, fraction) as floats.
    """
    match = re.search(r'ch([\d.]+)_cs([\d.]+)_fr([\d.]+)\.png', filename.replace("o","."))
    if match:
        charge = float(match.group(1).replace('o', '.'))
        conc = float(match.group(2).replace('o', '.'))
        fraction = float(match.group(3).replace('o', '.'))
        return charge, conc, fraction
    return None, None, None

def plot_images_on_grid(folder_path, fixed_axis, fixed_value):
    """
    Basic version: Plot images in a grid with two varying axes, one fixed axis.

    Parameters:
        folder_path (str): Path to the image folder.
        fixed_axis (str): One of 'charge', 'conc', or 'frac'.
        fixed_value (float): Fixed value for the chosen axis.
        extract_parameters (function): Function to extract (charge, conc, frac) from filename.
    """
    assert fixed_axis in ['charge', 'conc', 'frac'], "Invalid fixed axis"

    axis_map = {
        'charge': ('conc', 'frac', 'ch', fixed_value),
        'conc': ('charge', 'frac', 'cs', fixed_value),
        'frac': ('charge', 'conc', 'fr', fixed_value),
    }

    var1, var2, fixed_prefix, val = axis_map[fixed_axis]
    val_str = f"{fixed_prefix}{val:.2f}".replace(".", "o")
    files = [f for f in os.listdir(folder_path) if val_str in f]

    points = []
    for file in files:
        ch, cs, fr = extract_parameters(file)
        if ch is not None and cs is not None and fr is not None:
            values = {'charge': ch, 'conc': cs, 'frac': fr}
            points.append((values[var1], values[var2]))

    vals1 = sorted(set(p[0] for p in points))
    vals2 = sorted(set(p[1] for p in points))

    zoom_factor = min(1.0 / len(vals1), 1.0 / len(vals2)) * 0.65
    fig, ax = plt.subplots(figsize=(10, 8))

    ax.set_xticks(range(1, len(vals1) + 1))
    ax.set_yticks(range(1, len(vals2) + 1))
    ax.set_xlim([0.5, len(vals1) + 0.5])
    ax.set_ylim([0.5, len(vals2) + 0.5])

    for v1 in vals1:
        for v2 in vals2:
            params = {
                'charge': fixed_value if fixed_axis == 'charge' else v1 if var1 == 'charge' else v2,
                'conc': fixed_value if fixed_axis == 'conc' else v1 if var1 == 'conc' else v2,
                'frac': fixed_value if fixed_axis == 'frac' else v1 if var1 == 'frac' else v2,
            }
            filename = "ch{:.2f}_cs{:.2f}_fr{:.2f}.png".format(
                params['charge'], params['conc'], params['frac']
            ).replace(".", "o")
            filepath = os.path.join(folder_path, filename)
            if os.path.exists(filepath):
                img = plt.imread(filepath)
                image_box = OffsetImage(img, zoom=zoom_factor)
                ab = AnnotationBbox(image_box, (
                    vals1.index(v1) + 1, vals2.index(v2) + 1), frameon=False)
                ax.add_artist(ab)

    plt.tight_layout()
    return fig, ax  # You customize it further