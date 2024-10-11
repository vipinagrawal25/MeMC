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

file = sys.argv[1]
dirn = os.path.dirname(file)

nghst = 12
dtmp = lib.readHdf5(file, "pos")
Np = int(len(dtmp)/3)
pts_3d = dtmp.reshape(Np,3)
# lipA = np.loadtxt(dirn+"/lipA.txt")

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

outf=file.replace(".h5",".vtk")
tri_all = triangulate_solids(nbr_all)
lib.vtk_points(outf, pts_3d, tri_all)
lib.vtk_points_scalar(outf, pts_3d, lib.readHdf5(file,"lip"), name_scalar='lipid')