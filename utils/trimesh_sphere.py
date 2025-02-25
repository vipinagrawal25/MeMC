import trimesh
import matplotlib.pyplot as plt
from lib import *
import numpy as np
import sys

# Generate an icosphere with a specific subdivision level
mesh = trimesh.creation.icosphere(subdivisions=4, radius=20.0)
# print(mesh.faces)
# print(mesh.vertices)
Np = mesh.vertices.shape[0]
sort_tri = sort_simplices(mesh.faces)
cmlist, node_nbr = neighbours(Np, sort_tri)
node_nbr = sort_nbrs(mesh.vertices, Np, cmlist, node_nbr)
new_nbr = new_way_nbrs(Np, cmlist, node_nbr, nghst=12)
ncmlist = np.diff(cmlist)
write_hdf5(mesh.vertices, ncmlist, new_nbr, sys.argv[1])