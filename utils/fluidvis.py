import numpy as np
import h5py
import sys
import lib as lib
import os
import re

def make_triangles(start, nbrs, triangles):
    for i in range(len(nbrs)-1):
        triangles.append([start, nbrs[i], nbrs[i+1]])
    triangles.append([start,nbrs[0],nbrs[-1]])

def triangulate_solids(all_nbrs):
    tri_all = []
    for id_n,nbr_s in enumerate(all_nbrs):
        make_triangles(id_n, nbr_s, tri_all)
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

for file in sys.argv[1:]:
    # file = sys.argv[1]
    print(f"Processing {file}")
    outf=file.replace(".h5",".vtk")
    if os.path.exists(outf):
        print(f"{outf} already exists.")
        continue
    nghst = 12
    dtmp = lib.readHdf5(file, "pos")
    Np = int(len(dtmp)/3)
    pts_3d = dtmp.reshape(Np,3)

    Nnbr = lib.readHdf5(file, "node_nbr")
    Cmlst = lib.readHdf5(file, "cumu_list")

    all_pts =  np.full(Np, True, dtype=bool)
    all_id = all_pts

    nbr_all = []
    num_nbr_all = Cmlst[all_id]

    for k in range(Np):
        nnbr_id = num_nbr_all[k]
        st_idx = nghst*k
        nbr_all.append(Nnbr[st_idx:st_idx+nnbr_id])

    tri_all = triangulate_solids(nbr_all)
    lib.vtk_points(outf, pts_3d, tri_all)
    with h5py.File(file, 'a') as f:
        if 'cells' not in f:
            f.create_dataset('cells', data=tri_all)

    try:
        lip_data = lib.readHdf5(file, "lip")
        lib.vtk_points_scalar(outf, pts_3d, lip_data, name_scalar='lipid')
    except KeyError:
        print("Warning: 'lip' dataset not found in HDF5 file.")