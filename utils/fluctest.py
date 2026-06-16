import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import h5py
import glob

######################### FUNCTIONS ###############################
def h5_to_data(file,inputfile):
    with h5py.File(file,'r') as f:
        pos=f['pos'][()]
    pts=pos.reshape(-1,3)
    with h5py.File(inputfile,'r') as f:
        tri=f['triangles'][()]
    return pts, tri

def vertex_areas(pts,tri):
    N=pts.shape[0]
    areas=np.zeros(N)
    i=pts[tri[:,0]]
    j=pts[tri[:,1]]
    k=pts[tri[:,2]]
    cp=np.cross(j-i,k-i)
    tri_areas=0.5*np.linalg.norm(cp,axis=1)
    for l in range(3):
        np.add.at(areas,tri[:,l],tri_areas/3.0)
    return areas

inputfile=sys.argv[2]
dir=sys.argv[1]
LX=16*np.pi
LY=1*np.pi
Q_MAX=40

files=sorted(glob.glob(os.path.join(dir,"snap_*.h5")))
print(f"Found {len(files)} snaps.")
files = files[len(files)//3:]
print(f"Considering last {len(files)} snaps.")

spec={}
frames=0

def calc_hq2(pts,areas,LX,LY,Q_MAX):
    z=pts[:,2]
    h=z-np.mean(z)
    x,y=pts[:,0],pts[:,1]
    qs={}
    for nx in range(-Q_MAX,Q_MAX+1):
        for ny in range(-Q_MAX,Q_MAX+1):
            if nx==0 and ny==0:
                continue
            qx,qy=2*np.pi*nx/LX,2*np.pi*ny/LY
            q=np.sqrt(qx**2+qy**2)
            phase=qx*pts[:,0]+qy*pts[:,1]
            hq=np.sum(h*np.exp(-1j*phase)*areas)
            hq2=np.abs(hq)**2/(LX*LY)
            qround=np.round(q,3)
            qs.setdefault(qround,[]).append(hq2)
    return {q: np.mean(vals) for q,vals in qs.items()}

for file in files:
    pts,tri=h5_to_data(file,inputfile)
    areas=vertex_areas(pts,tri)
    framespec=calc_hq2(pts,areas,LX,LY,Q_MAX)
    for q,val in framespec.items():
        spec[q]=spec.get(q,0)+val
    frames+=1
    if frames%10==0:
        print(f"At {frames}/{len(files)} frames.")


qvals=sorted(spec.keys())
hqmean=[spec[q]/frames for q in qvals]

kBT=1.0
kappa=4.0
qvals=np.asarray(qvals)
q0=hqmean[0]/qvals[0]**(-2)
plt.loglog(qvals,hqmean,'o-',label="data")
plt.loglog(qvals,q0*qvals**-2,'--',label="q^-2")
plt.legend()
plt.savefig("fluctuation_spectrum.png")
