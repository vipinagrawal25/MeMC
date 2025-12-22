import numpy as np
import h5py
import sys
import lib as lib
import os
import re
import argparse

def check_hdf5_structure(filename):
    """Check and print the structure of an HDF5 file"""
    print(f"HDF5 file structure for {filename}:")
    try:
        with h5py.File(filename, 'r') as f:
            def print_structure(name, obj):
                if isinstance(obj, h5py.Dataset):
                    print(f"  Dataset: {name} - shape: {obj.shape}, dtype: {obj.dtype}")
                elif isinstance(obj, h5py.Group):
                    print(f"  Group: {name}")
            f.visititems(print_structure)
    except Exception as e:
        print(f"Error reading HDF5 file: {e}")
    print()

def make_triangles(start, nbrs, triangles, points=None, box_size=None):
    for i in range(len(nbrs)-1):
        if points is not None and box_size is not None:
            # Check if the triangle would cross periodic boundaries
            p1, p2, p3 = points[start], points[nbrs[i]], points[nbrs[i+1]]
            if any(abs(p2-p1) > box_size/2) or any(abs(p3-p1) > box_size/2) or any(abs(p3-p2) > box_size/2):
                continue
        triangles.append([start, nbrs[i], nbrs[i+1]])
    
    # Check the last triangle (connecting first and last neighbors)
    if points is not None and box_size is not None:
        p1, p2, p3 = points[start], points[nbrs[0]], points[nbrs[-1]]
        if not any(abs(p2-p1) > box_size/2) and not any(abs(p3-p1) > box_size/2) and not any(abs(p3-p2) > box_size/2):
            triangles.append([start, nbrs[0], nbrs[-1]])

def triangulate_solids(all_nbrs, points=None, box_size=None):
    tri_all = []
    for id_n,nbr_s in enumerate(all_nbrs):
        make_triangles(id_n, nbr_s, tri_all, points, box_size)
    return tri_all

def extract_parameters(filename):
    # filename = filename.split("/snap")[0].replace("/","_")
    if "snap" in filename:
        match = re.search(r'ch([\d.]+)/cs([\d.]+)/fr([\d.]+)/', filename.replace("o","."))
    else:
       match = re.search(r'ch([\d.]+)_cs([\d.]+)_fr([\d.]+)\.h5', filename.replace("o","."))
    if match:
        charge = float(match.group(1).replace('o', '.'))
        conc = float(match.group(2).replace('o', '.'))
        fraction = float(match.group(3).replace('o', '.'))
        return charge, conc, fraction
    return None, None, None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert HDF5 files to VTK format')
    parser.add_argument('files', nargs='+', help='HDF5 files to process')
    parser.add_argument('--inspect-only', action='store_true', 
                       help='Only inspect HDF5 structure without processing')
    
    args = parser.parse_args()
    
    for file in args.files:
        print(f"Processing {file}")
        check_hdf5_structure(file)
        
        if args.inspect_only:
            continue
            
        outf=file.replace(".h5",".vtk")
        if os.path.exists(outf):
            print(f"{outf} already exists.")
            continue
        nghst = 12
        
        try:
            dtmp = lib.readHdf5(file, "pos")
        except KeyError:
            print("Error: 'pos' dataset not found in HDF5 file. Cannot process without position data.")
            continue
            
        Np = int(len(dtmp)/3)
        pts_3d = dtmp.reshape(Np,3)

        try:
            Nnbr = lib.readHdf5(file, "node_nbr")
            Cmlst = lib.readHdf5(file, "cumu_list")
        except KeyError as e:
            print(f"Warning: Required dataset not found in HDF5 file: {e}")
            print("Available datasets in file:")
            with h5py.File(file, 'r') as f:
                def print_datasets(name, obj):
                    if isinstance(obj, h5py.Dataset):
                        print(f"  - {name}: shape {obj.shape}, dtype {obj.dtype}")
                f.visititems(print_datasets)
            print("Skipping file due to missing required datasets.")
            continue

        all_pts =  np.full(Np, True, dtype=bool)
        all_id = all_pts

        nbr_all = []
        num_nbr_all = Cmlst[all_id]

        for k in range(Np):
            nnbr_id = num_nbr_all[k]
            st_idx = nghst*k
            nbr_all.append(Nnbr[st_idx:st_idx+nnbr_id])

        # Try to read box size from the HDF5 file
        # Insert this function (e.g., after extract_parameters)
        def estimate_box_size_from_positions(pts):
            """Estimate cubic box size from particle positions. Returns box_size or None."""
            try:
                mins = pts.min(axis=0)
                maxs = pts.max(axis=0)
                ranges = maxs - mins
                box_size = float(np.max(ranges))
                if box_size == 0 or not np.isfinite(box_size):
                    print("Warning: unable to estimate box size from positions.")
                    return None
                # tiny tolerance to avoid strict half-box checks failing due to numerical precision
                box_size *= 1.0
                print(f"Estimated box size from positions: {box_size}")
                return box_size
            except Exception as e:
                print(f"Warning: failed to estimate box size automatically: {e}")
                return None

                # Replace the original placeholder code with this call:
                # Estimate box size automatically from particle positions
        box_size = estimate_box_size_from_positions(pts_3d)
            
        tri_all = triangulate_solids(nbr_all, pts_3d if box_size is not None else None, box_size)
        lib.vtk_points(outf, pts_3d, tri_all)
        with h5py.File(file, 'a') as f:
            if 'cells' not in f:
                f.create_dataset('cells', data=tri_all)

        try:
            lip_data = lib.readHdf5(file, "lip")
            lib.vtk_points_scalar(outf, pts_3d, lip_data, name_scalar='lipid')
        except KeyError:
            print("Warning: 'lip' dataset not found in HDF5 file.")
            frac = None
            try:
                _, _, frac = extract_parameters(file)
            except Exception as e:
                print(f"Could not extract fraction parameter: {e}")
            if frac == 0:
                print("Creating lipid data for frac=0 (all zeros)")
                lip_data = np.zeros(Np)
            elif frac == 1:
                print("Creating lipid data for frac=1 (all ones)")
                lip_data = np.ones(Np)
            else:
                lip_data = None
            if lip_data is not None:
                lib.vtk_points_scalar(outf, pts_3d, lip_data, name_scalar='lipid')
        
        print(f"Successfully processed {file} -> {outf}")
        print()