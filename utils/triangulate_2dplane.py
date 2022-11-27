#+begin_src python :session py2

import numpy as np
import matplotlib.pyplot as plt
import h5py, sys
from scipy.spatial import Delaunay

#+end_src

#+RESULTS:



#+begin_src python :session py2

def read_data(filename):
    # pos = h5py.File(filename)["pos"][()]
    # pos = np.asarray(pos)
    # Np = int(len(pos)/2)
    # pts = pos.reshape(Np,2)
    pts = np.loadtxt(filename)
    Np = len(pts);
    return Np, pts

def sort_simplices(cells):
    lsimples = len(cells)
    nsimplices = np.asarray([], dtype=np.int32)
    for scles in cells:
        nscles = np.sort(scles)
        nsimplices = np.hstack([nsimplices, nscles])
        nsimplices = np.hstack([nsimplices, [nscles[1], nscles[2], nscles[0]]])
        nsimplices = np.hstack([nsimplices, [nscles[2], nscles[0], nscles[1]]])
        nsimplices = np.hstack([nsimplices, [nscles[0], nscles[2], nscles[1]]])
        nsimplices = np.hstack([nsimplices, [nscles[1], nscles[0], nscles[2]]])
        nsimplices = np.hstack([nsimplices, [nscles[2], nscles[1], nscles[0]]])
    nsimplices = nsimplices.reshape(lsimples*6, 3)
    nsimplices = np.asarray(sorted(nsimplices, key=lambda x: (x[0], x[1])))
    return nsimplices

def neighbours(Np, simpl):
    r1=simpl[:,0]
    r2=simpl[:,1]
    r3=simpl[:,2]
    lst=np.zeros(Np,dtype=int)
    cumlst=np.zeros(Np+1,dtype=int)
    for i in range(0, Np):
        lst[i]=len(r1[r1==i])/2

    cumlst[1:] = np.cumsum(lst)
    node_neighbour = np.zeros(cumlst[-1],dtype=int)
    for i in range(0, cumlst[-1], 1):
        node_neighbour[i]=r2[2*i]
    return cumlst,node_neighbour



def write_hdf5(R, cmlst, node_nbr,  posfile, file):
    if file.split(".")[-1]=="h5":
        pass
    else:
        file=file+".h5"
    hf = h5py.File(posfile,'w')
    hf.create_dataset('pos',data=R.reshape(-1))
    hf.close()

    hf = h5py.File(file,'w')
    hf.create_dataset('cumu_list',data=cmlst.astype(np.int32))
    hf.create_dataset('node_nbr',data=node_nbr.astype(np.int32))
    hf.close()

# def write_file(pts_cart, cmlist, node_nbr):
#     file = open("../Examples/pts.bin", "wb")
#     file.write(pts_cart)
#     file.close()

#     file = open("../Examples/mesh.bin", "wb")
#     file.write(cmlist)
#     file.write(node_nbr)
#     file.close()

#+end_src

#+RESULTS:

#+begin_src python :session py2 :results output

file = 'points__.dat'
def read_simplices(file):
    simplices = []
    data = np.loadtxt(file, dtype=int)
    for i in range(0,len(data)):
        simplices.append([data[i,0], data[i,1], data[i,2]])
    return simplices

def new_way_nbrs(cmlist, node_nbr, nghst=12):
    new_nbr = np.zeros(nghst*Np, dtype=int)
    new_nbr[:] = -1
    print(new_nbr)
    for ip in range(0, Np):
        nbrs = node_nbr[cmlist[ip]:cmlist[ip+1]]
        num_nbr = -(cmlist[ip]-cmlist[ip+1])
        angles = []
        for i in nbrs:
            dx = pts[i,0]-pts[ip,0]
            if(dx >= 0.5): dx = dx - 1.0
            if(dx < -0.5): dx = dx + 1.0
            dy = pts[i,1]-pts[ip,1]
            if(dy >= 0.5): dy = dy - 1.0
            if(dy < -0.5): dy = dy + 1.0
            angles.append(np.arctan2(dy, dx))
        angles = np.asarray(angles)
        idx = np.where(angles < 0)
        angles[idx] = angles[idx] + 2*np.pi
        sort = np.argsort(angles)
        nnbrs = nbrs[sort]
        st_idx = int(ip*nghst); end_idx = int(ip*nghst + num_nbr)
        new_nbr[st_idx:end_idx] = nnbrs[:]

    return new_nbr



Np, pts = read_data(file)
pts_3d = np.zeros(shape=(Np, 3), dtype=float)
pts_3d[:,0] = pts[:,0]
pts_3d[:,1] = pts[:,1]
# tri = Delaunay(pts, furthest_site = False)
simplices = read_simplices('simplices__.dat')
# ng = 128
sort_tri = sort_simplices(simplices)
cmlist, node_nbr = neighbours(Np, sort_tri)
new_nbr = new_way_nbrs(cmlist, node_nbr, nghst=12)
ncmlist = np.diff(cmlist)
write_hdf5(pts_3d,  ncmlist, new_nbr, 'conf/inp_pos.h5', 'conf/nbrs.h5')

#+end_src

#+RESULTS:
: [-1 -1 -1 ... -1 -1 -1]

#+begin_src python :session py2 :results output

fig, ax = plt.subplots()
# ax.plot(pts[:,0], pts[:,1], '.')
i = 32
print(ncmlist[i])
x = pts[i,0]
if (x < 0.5): x = 1 + x
ax.plot(x, pts[i,1], 'o')
for j in new_nbr[i*12:i*12 + ncmlist[i]]:
    x = pts[j,0]
    if (x < 0.5): x = 1 + x
    ax.plot(x, pts[j,1], 's')

fig, ax = plt.subplots()
for i in range(0, 6):
    ax.plot(i, 0, '-o')

plt.show()

#+end_src

#+RESULTS:
: 6
