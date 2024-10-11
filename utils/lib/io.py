import h5py as h5py
import numpy as np

def writeHdf5(R, cmlst, node_nbr,  posfile, file):
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


def readHdf5(fn, grpname):
    hf = h5py.File(fn,'r')
    data = hf[grpname][()]
    hf.close()
    return data