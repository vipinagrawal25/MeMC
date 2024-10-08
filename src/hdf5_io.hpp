#ifndef HDF5_IO_HPP
#define HDF5_IO_HPP
#include <hdf5.h>
#include <unistd.h>
#include "vector.hpp"

void hdf5_io_read_double(double *Pos, string input_file, string);
void hdf5_io_read_mesh(int *cmlist, int *node_nbr, string input_file);
void hdf5_io_write_mesh(int *cmlist,
        int *node_nbr, int N, int ng, string output_file);
void io_dump_config(double *Pos, int N, char *);
void io_read_config(double *Pos, int N, char *);
void io_dump_config_ascii(double *Pos, int N, char *);
int hdf5_io_get_Np(string input_file, string dset_name);
void hdf5_io_delete(string filename);
//----//
template <typename T>
void hdf5_io_write(T *Pos, int N, string input_file, string dset_name){
    ///  @brief hdf5 io for  the position 
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param N number of points 
    ///  @param input_file File name to dump all the co-ordinate
    hid_t   file_id, dset1, space_id;  /* identifiers */
    herr_t  status;
    hsize_t  dims; 

    hid_t h5_type;
    if (typeid(T) == typeid(double)){
        h5_type = H5T_NATIVE_DOUBLE;
    } else if (typeid(T) == typeid(int)){
        h5_type = H5T_NATIVE_INT;
    } else if (typeid(T) == typeid(bool)){
        h5_type=H5T_NATIVE_HBOOL;
    }else{
        cerr << "Unsupported data type!" << endl;
        return;
    }

    if(access(input_file.c_str(),F_OK)!=0){
        file_id = H5Fcreate (input_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }else{
        /* Open an existing file. */
        file_id = H5Fopen (input_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }

    // Check if the dataset already exists
    if (H5Lexists(file_id, dset_name.c_str(), H5P_DEFAULT) > 0) {
        H5Ldelete(file_id, dset_name.c_str(), H5P_DEFAULT);
    }

    dims = N;
    
    space_id = H5Screate_simple (1, &dims, NULL);
    dset1 = H5Dcreate2(file_id, dset_name.c_str(), h5_type, space_id, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset1, h5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                Pos);
    status = H5Dclose (dset1);
    status = H5Sclose (space_id);
    status = H5Fclose (file_id);
    
    if(status != 0){
      cerr << "File close failed" << endl;
    }
}
#endif