import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os
import glob
import sys

def plot_and_save_vtk_files(folder_path, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    
    custom_cmap = LinearSegmentedColormap.from_list("blue_white", [(1, 1, 1), (0, 0, 1)], N=256)
    
    vtk_files = glob.glob(os.path.join(folder_path, "*.vtk"))
    
    if not vtk_files:
        print("No .vtk files found in the folder.")
        return

    for vtk_file in vtk_files:
        base_name = os.path.basename(vtk_file).replace(".vtk", ".png")
        output_file = os.path.join(output_folder, base_name)

        # Check if the output file already exists
        if os.path.exists(output_file):
            print(f"Skipping: {vtk_file} (Output file already exists)")
            continue

        print(f"Processing: {vtk_file}")

        # Load the .vtk file
        mesh = pv.read(vtk_file)

        scalars = None
        if mesh.point_data:
            scalars = list(mesh.point_data.keys())[0]
        elif mesh.cell_data:
            scalars = list(mesh.cell_data.keys())[0]

        # Set up the plotter
        plotter = pv.Plotter(off_screen=True)  # Use off_screen for saving images without showing
        plotter.add_mesh(mesh, cmap=custom_cmap, show_edges=False)

        # Remove scalar bar if it exists
        if plotter.scalar_bars:
            plotter.remove_scalar_bar()

        plotter.screenshot(output_file, window_size=[800, 800], scale=1.0)
        plotter.close()
        # Show the plotter window (for debugging purposes)
        print(f"Saved: {output_file}")
    
# Example usage
folder_path = sys.argv[1]  # Replace with the path to your folder containing .vtk files
output_folder = sys.argv[1]  # Replace with the ipath to save the images
plot_and_save_vtk_files(folder_path, output_folder)