import h5py as h5py
import numpy as np

def write_hdf5(R, cmlst, node_nbr,  file):
    if file.split(".")[-1]=="h5":
        pass
    else:
        file=file+".h5"
    hf = h5py.File(file,'w')
    hf.create_dataset('pos',data=R.reshape(-1))
    hf.create_dataset('cumu_list',data=cmlst.astype(np.int32))
    hf.create_dataset('node_nbr',data=node_nbr.astype(np.int32))
    hf.close()

def readHdf5(fn, grpname):
    hf = h5py.File(fn,'r')
    data = hf[grpname][()]
    hf.close()
    return data

def save_cell_dataset(dataset, filename):
    if not filename.endswith('.h5'):
        filename = filename + '.h5'
    
    with h5py.File(filename, 'w') as hf:
        # for key, data in dataset.items():
        hf.create_dataset("cells", data=dataset)

    print(f"Cell dataset saved to {filename}")

def write_pos_cells(R, cells, file):
    if file.split(".")[-1]=="h5":
        pass
    else:
        file=file+".h5"
    hf = h5py.File(file,'w')
    hf.create_dataset('pos',data=R.reshape(-1))
    hf.create_dataset('cells',data=cells)
    hf.close()