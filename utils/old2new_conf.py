import numpy as np
import h5py, sys
#
def write_hdf5(R, cmlst, node_nbr,  posfile, file):
 if file.split(".")[-1]=="h5":
     pass
 else:
     file=file+".h5"
 hf = h5py.File(posfile,'w')
 hf.create_dataset('pos',data=R.reshape(-1))
 # hf.close()
 # hf = h5py.File(file,'w')
 hf.create_dataset('cumu_list',data=cmlst.astype(np.int32))
 hf.create_dataset('node_nbr',data=node_nbr.astype(np.int32))
 hf.close()

def read_hdf5(posfile, file):
 if file.split(".")[-1]=="h5":
     pass
 else:
     file=file+".h5"
 hf = h5py.File(posfile,'r')
 R = hf["pos"][()]
 Np = int(len(R)/3)
 hf.close()

 hf = h5py.File(file,'r')
 cmlst = hf["cumu_list"][()]
 node_nbr = hf["node_nbr"][()]
 hf.close()
 return R, cmlst, node_nbr, Np

def new_way_nbrs(cmlist, node_nbr, nghst=12):
 new_nbr = np.zeros(nghst*Np, dtype=int)
 new_nbr[:] = -1
 print(new_nbr)
 for ip in range(0, Np):
     nbrs = node_nbr[cmlist[ip]:cmlist[ip+1]]
     num_nbr = -(cmlist[ip]-cmlist[ip+1])
     nnbrs = nbrs
     st_idx = int(ip*nghst); 
     end_idx = int(ip*nghst + num_nbr)
     new_nbr[st_idx:end_idx] = nnbrs[:]
 return new_nbr
pos_file = sys.argv[1]
mesh_file = sys.argv[2]
out_file = sys.argv[3]
pts_3d, cmlist, node_nbr, Np = read_hdf5(pos_file, mesh_file)
print(Np)
new_nbr = new_way_nbrs(cmlist, node_nbr, nghst=8)
ncmlist = np.diff(cmlist)
write_hdf5(pts_3d,  ncmlist, new_nbr, out_file, out_file)