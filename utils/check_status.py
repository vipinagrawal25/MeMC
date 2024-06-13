import numpy as np
import matplotlib.pyplot as plt
import h5py, sys

def read_data(filename, which_file):
    pos = h5py.File(filename)["pos"][()]
    pos = np.asarray(pos)
    if(which_file == 'start'):
        Np = int(len(pos)/2)
        pts_sph = pos.reshape(Np,2)
        pts_cart = np.asarray([[np.sin(theta)*np.cos(phi), 
            np.sin(theta)*np.sin(phi), 
            np.cos(theta)] for theta, phi in pts_sph]
            ) 
    elif(which_file == 'memc'):
        Np = int(len(pos)/3)
        pts_cart = pos.reshape(Np,3)
        pts_sph = np.asarray([[np.sqrt(x**2 + y**2 + z**2), 
            np.arctan2(np.sqrt(x**2 + y**2), z), 
            np.arctan2(y, x)] for x, y, z in pts_cart]
            ) 
    return Np, pts_sph, pts_cart
    
which_file = sys.argv[1]
file = sys.argv[2]
fig, ax = plt.subplots()
Np, pts_sph, pts_cart = read_data(file, which_file)
np.savetxt('test.dat', pts_sph/(2*np.pi), fmt='%.16f')
if(which_file == 'start'):
    ax.plot((pts_sph[0:64,0]), pts_sph[0:64,1], 'o')
    ax.plot((pts_sph[64:,0]), pts_sph[64:,1], 'o')
if(which_file == 'memc'):
    ax.plot(np.cos(pts_sph[:,1]), pts_sph[:,2], 's')

# ax.set(xlim = [-1.1,1.1], ylim=[-0.1,6.36], 
#         xlabel='cos(theta)', ylabel='phi')

plt.show()
