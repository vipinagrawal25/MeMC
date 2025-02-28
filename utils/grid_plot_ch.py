import os
import re
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
home_dir = os.path.expanduser("~")
plt.style.use(os.path.join(home_dir, '.matplotlibrc'))
import numpy as np
import sys

def extract_parameters(filename):
    match = re.search(r'ch([\d.]+)_cs([\d.]+)_fr([\d.]+)\.png', filename.replace("o","."))
    if match:
        charge = float(match.group(1).replace('o', '.'))
        conc = float(match.group(2).replace('o', '.'))
        fraction = float(match.group(3).replace('o', '.'))
        return charge, conc, fraction
    return None, None, None

def plot_images_on_grid(folder_path, charge=0.03, concs=None, fracs=None):
    # Get all .png files in the folder
    charge_str = "ch{:.2f}".format(charge).replace(".", "o")
    files = [f for f in os.listdir(folder_path) if charge_str in f]
    # Extract unique charge and concentration values
    points = []
    for file in files:
        _, conc, frac = extract_parameters(file)
        if frac is not None and conc is not None:
            points.append((conc, frac))

    if concs is None:
        concs = sorted(set(pt[0] for pt in points))
    if fracs is None:
        fracs = sorted(set(pt[1] for pt in points))

    print(concs, fracs)
    # Calculate zoom factor based on the number of images
    zoom_factor = min(1.0 / len(concs), 1.0 / len(fracs))*0.7

    # print(charges, concs)
    # Create a grid for plotting
    fig, ax = plt.subplots(figsize=(10, 8))

    # Adjust subplot parameters to reduce space between axes and plots
    
    # Remove the top and right lines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Remove all axis elements
    ax.set_xticks(range(1, len(concs) + 1))
    ax.set_xlim([0.5,len(concs)+0.5])
    ax.set_ylim([0.3,len(fracs)+0.5])

    ax.set_xticklabels([f"{c:.2f}" for c in concs])
    ax.set_yticks(range(1, len(fracs) + 1))
    ax.set_yticklabels([f"{c:.2f}" for c in fracs])
    ax.set_xlabel("Salt concentration (mM)")
    ax.set_ylabel("Fraction")

    # Add xticklabels on top with values = radius*sqrt(c)/0.304
    radius = 11.3  # Define the radius value
    top_xticklabels = [f"{radius * np.sqrt(c) / 0.304:.2f}" for c in concs]
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    ax_top.set_xticks(range(1, len(concs) + 1))
    ax_top.set_xticklabels(top_xticklabels)
    ax_top.set_xlabel(r"${R}/{\lambda_D}$")

    # Remove short ticks
    ax.tick_params(axis='both', which='minor', length=0)
    # Create the output directory if it does not exist
    output_dir = folder_path.replace("/latest_files", "/plots/")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # Plot images at appropriate locations
    for frac in fracs:
        for conc in concs:
            file = "ch{:.2f}_cs{:.2f}_fr{:.2f}".format(charge, conc, frac)
            filepath = os.path.join(folder_path, file.replace(".", "o") + ".png")
            # Load the image
            if os.path.exists(filepath):
                img = plt.imread(filepath)
                # Adjust zoom to control image size
                image_box = OffsetImage(img, zoom=zoom_factor)
                bboxprops = dict(edgecolor='black', linewidth=1)  # Add line around the image
                ab = AnnotationBbox(image_box, (concs.index(conc) + 1, fracs.index(frac) + 1), frameon=False, bboxprops=bboxprops)
                ax.add_artist(ab)
    output_filename = output_dir+f"ch{charge:.2f}".replace(".", "o") + ".png"
    title = f"ch={charge:.2f}"

    plt.savefig(output_filename, bbox_inches='tight')
    
# Example usage
folder_path = sys.argv[1]+"/latest_files"  # Replace with your folder path
charges = sorted(set(extract_parameters(f)[0] for f in os.listdir(folder_path) if extract_parameters(f)[0] is not None))
print(charges)
#
for charge in charges:
    plot_images_on_grid(folder_path,charge)