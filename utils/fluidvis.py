import numpy as np
import h5py
import sys
import lib as lib
import os

def make_triangles(start, nbrs, triangles):
    for i in range(len(nbrs)-1):
        triangles.append([start, nbrs[i], nbrs[i+1]])
    triangles.append([start,nbrs[0],nbrs[-1]])

def triangulate_solids(all_nbrs):
    tri_all = []
    for id_n,nbr_s in enumerate(all_nbrs):
        make_triangles(id_n, nbr_s, tri_all)
    return tri_all

def estimate_box_size(pts_3d):
    """Estimate box size from point coordinates"""
    ranges = np.ptp(pts_3d, axis=0)  # peak-to-peak (max - min) for each dimension
    return np.max(ranges) * 1.1  # Add 10% margin

def filter_wrapped_neighbors(pts_3d, vertex_idx, neighbor_indices, box_size):
    """Filter out neighbors that are wrapped around periodic boundaries"""
    if box_size <= 0:
        return neighbor_indices  # No filtering if box size unknown
    
    vertex_pos = pts_3d[vertex_idx]
    filtered_neighbors = []
    
    for nbr_idx in neighbor_indices:
        nbr_pos = pts_3d[nbr_idx]
        distance = np.linalg.norm(vertex_pos - nbr_pos)
        
        # If distance is greater than half box size, likely wrapped
        if distance <= box_size / 2:
            filtered_neighbors.append(nbr_idx)
        else:
            print(f"Filtering wrapped neighbor: vertex {vertex_idx} -> neighbor {nbr_idx} (distance: {distance:.3f})")
    
    return filtered_neighbors

for file in sys.argv[1:]:
    print(f"Processing {file}")
    dirn = os.path.dirname(file)
    outf=file.replace(".h5",".vtk")
    if os.path.exists(outf):
        print(f"{outf} already exists.")
        continue
    nghst = 12
    dtmp = lib.readHdf5(file, "pos")
    Np = int(len(dtmp)/3)
    pts_3d = dtmp.reshape(Np,3)

    # Estimate box size for periodic boundary detection
    box_size = estimate_box_size(pts_3d)
    print(f"Estimated box size: {box_size:.3f}")

    Nnbr = lib.readHdf5(file, "node_nbr")
    Cmlst = lib.readHdf5(file, "cumu_list")

    all_pts =  np.full(Np, True, dtype=bool)
    all_id = all_pts

    nbr_all = []
    num_nbr_all = Cmlst[all_id]

    for k in range(Np):
        nnbr_id = num_nbr_all[k]
        st_idx = nghst*k
        raw_neighbors = Nnbr[st_idx:st_idx+nnbr_id]
        # Filter out wrapped neighbors
        filtered_neighbors = filter_wrapped_neighbors(pts_3d, k, raw_neighbors, box_size)
        nbr_all.append(filtered_neighbors)

    tri_all = triangulate_solids(nbr_all)
    lib.vtk_points(outf, pts_3d, tri_all)