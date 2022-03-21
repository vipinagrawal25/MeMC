import numpy as np
import h5py
import vtkio as vtk
import sys

inf = sys.argv[1]
conf_f = sys.argv[2]
outf = sys.argv[3]





def read_data(filename):
    f = h5py.File(filename, 'r')
    pos = f["pos"][()]
    dim = int(len(pos)/3)
    pos = pos.reshape(dim,3)
    return pos
 
triangles = h5py.File(conf_f)["triangles"][()]

pos = read_data(inf)

vtk.vtk_points(outf, pos, triangles)


