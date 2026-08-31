import numpy as np
import numpy.linalg as la
import sys
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import h5py

def read_data(filename,inputfile):
    pos=h5py.File(filename)["pos"][()]
    pos=np.asarray(pos)
    Np=int(len(pos)/3) #3D position coordinates
    pts=pos.reshape(Np,3) #we get the cartesian x,y coords
    pos=pos.reshape(-1,3)
    tri=h5py.File(inputfile)["triangles"][()] #triangulation data
    x=pts[:,0]
    y=pts[:,1]
    z=pts[:,2] #this is what we will be needing for the height field data
    mean_ht=np.mean(z)
    h=z-mean_ht
    return Np,pos,x,y,h,tri #we get the height field, x-y positions, and triangulation data

def get_ht_field(pts):
    z_cart=pts_cart[:,2]
    return z_cart

#call this program as "fluctuation_spectra.py [snap.h5] [input.h5]"
file=sys.argv[1]
inputfile=sys.argv[2]
Np,pos,x,y,h,tri=read_data(file,inputfile)
np.savetxt('heights.dat',np.column_stack([x,y,h]))
np.savetxt('triangles.dat',tri)

#check once by plotting
fig,ax=plt.subplots()
sc=ax.scatter(x,y,c=h,s=10)
plt.colorbar(sc,ax=ax,label="h(x,y)")
plt.savefig("height_field.png")
plt.close()

#now we need to compute the Fourier spectra from the height field
def find_vertex_areas(pos,triangles):
    N=len(pos)
    areas=np.zeros(N)
    for tri in triangles:
        i,j,k=tri #each triangle
        e1=pos[j]-pos[i]
        e2=pos[k]-pos[i]
        area=0.5*la.norm(np.cross(e1,e2))
        areas[i]+=area/3.0
        areas[j]+=area/3.0
        areas[k]+=area/3.0
    return areas

def hq(x,y,h,areas,qx,qy):
    exp=np.exp(1j*(qx*x+qy*y))
    return np.sum(phase*h*areas)

areas=find_vertex_areas(pos,tri)
total_area=areas.sum()
x0=pos[:,0]
y0=pos[:,1]
Lx=x0.max()-x0.min()
Ly=y0.max()-y0.min()
print(f"Lx={Lx:.4f}  Ly={Ly:.4f}")
print(f"Total area={total_area:.4f}  expected={Lx*Ly:.4f}")
