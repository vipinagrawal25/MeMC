/* #include <hdf5/serial/hdf5.h> */
// #include "global.h"
#include "hdf5_io.hpp"
#include <cstdio>
#include "misc.hpp"
 /**  
 *  @brief hdf5 IO for the mesh  
 *  
 */

int hdf5_io_get_Np(string input_file, string dset_name){
   hid_t file_id,dataset_id;  /* identifiers */
   herr_t status;
   // string dset_name="pos";

   if(access(input_file.c_str(),F_OK)!=0){
      fprintf(stderr, "The configuration file does not exit\n");
      exit(1);
   }

   file_id = H5Fopen(input_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   dataset_id = H5Dopen(file_id, dset_name.c_str(), H5P_DEFAULT);

   hid_t dspace = H5Dget_space(dataset_id);
   const int ndims = H5Sget_simple_extent_ndims(dspace);
   hsize_t dims[2];  // Support up to 2D arrays
   H5Sget_simple_extent_dims(dspace, dims, NULL);

   int total_points;
   if (ndims == 1) {
       total_points = dims[0];
   } else if (ndims == 2) {
       total_points = dims[0] * dims[1];  // For [N,3] this will be 3N
   } else {
       fprintf(stderr, "Unsupported number of dimensions: %d\n", ndims);
       H5Sclose(dspace);
       H5Dclose(dataset_id);
       H5Fclose(file_id);
       exit(1);
   }

   H5Sclose(dspace);
   H5Dclose(dataset_id);
   H5Fclose(file_id);

   return total_points;
}

void hdf5_io_read_double(double *Data, string input_file, string dset_name){
    ///  @brief Read from the hdf5 file
    ///  @param Data array containing co-ordinates of all the particles
    ///  @param input_file File name from which co-ordinate will be read
    ///  @note Handles both 1D arrays and 2D arrays with second dimension of 3

    hid_t file_id, dataset_id, dataspace;  /* identifiers */
    herr_t status;
    int ndims;
    hsize_t dims[2];  // Array to hold dimensions

    if(access(input_file.c_str(),F_OK)!=0){
        std::cerr << "Error: The configuration file does not exist\n";
        std::exit(EXIT_FAILURE);
    }
    
    /* Open an existing file. */
    file_id = H5Fopen(input_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        cerr << "Error: Could not open file " << input_file << "\n";
        exit(EXIT_FAILURE);
    }

    dataset_id = H5Dopen(file_id, dset_name.c_str(), H5P_DEFAULT);
    if (dataset_id < 0) {
        std::cerr << "Error: Could not open dataset " << dset_name << " in file " << input_file << "\n";
        H5Fclose(file_id);
        std::exit(EXIT_FAILURE);
    }

    // Get the dataspace and check dimensions
    dataspace = H5Dget_space(dataset_id);
    ndims = H5Sget_simple_extent_ndims(dataspace);
    status = H5Sget_simple_extent_dims(dataspace, dims, NULL);
    
    if (ndims == 2) {
        // For 2D dataset, verify second dimension is 3
        if (dims[1] != 3) {
            std::cerr << "Error: For 2D dataset, second dimension must be 3 for (x,y,z), got " 
                      << dims[1] << "\n";
            H5Sclose(dataspace);
            H5Dclose(dataset_id);
            H5Fclose(file_id);
            std::exit(EXIT_FAILURE);
        }
    }

    // Read the data
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Data);
    if (status < 0) {
        std::cerr << "Error: Failed to read dataset " << dset_name << "\n";
        H5Sclose(dataspace);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
        std::exit(EXIT_FAILURE);
    }

    H5Sclose(dataspace);
    status = H5Dclose(dataset_id);
    if (status < 0) {
        std::cerr << "Error: Failed to close dataset " << dset_name << "\n";
    }
    status = H5Fclose(file_id);
    if (status < 0) {
        std::cerr << "Error: Failed to close file " << input_file << "\n";
    }
}




void hdf5_io_write_mesh(int *cmlist, int *node_nbr, int N, int ng, string output_file){

    ///  @brief Read the mesh from the hdf5 file
    ///  @param cmlist array containing the number of neighbours for each particle  
    ///  @param node_nbr array containing the list of neighbours for each particle  
    ///  @param input_file File name from which co-ordinate will be read
    /// 

    hid_t   file_id, dset1, dataset_id, space_id;  /* identifiers */
    herr_t  status;
    int size_mesh; 
    hsize_t          dims; 

    size_mesh = ng*N;

    if(access(output_file.c_str(),F_OK)!=0){
        file_id = H5Fcreate (output_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }else{

        file_id = H5Fopen (output_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }


    dims = N;
    space_id = H5Screate_simple (1, &dims, NULL);
    // cout << space_id << endl;
    dset1 = H5Dcreate2(file_id, "/cumu_list", H5T_NATIVE_INT, space_id, H5P_DEFAULT,
            H5P_DEFAULT, H5P_DEFAULT);
    // cout << dset1 << endl;
    if (dset1 < 0) {
        fprintf(stderr, "Failed to create dataset for cumu_list.\n");
        H5Sclose(space_id);
        H5Fclose(file_id);
    }

    status = H5Dwrite (dset1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            cmlist);
    status = H5Dclose (dset1);
    status = H5Sclose (space_id);

    dims = size_mesh;
    space_id = H5Screate_simple (1, &dims, NULL);
    dset1 = H5Dcreate2(file_id, "/node_nbr", H5T_NATIVE_INT, space_id, H5P_DEFAULT,
            H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            node_nbr);
    status = H5Dclose (dset1);
    status = H5Sclose (space_id);
    status = H5Fclose(file_id);
    if(status != 0){
        fprintf(stderr, "file close failed\n");
    }
}


void hdf5_io_read_mesh(int *cmlist, int *node_nbr,  string input_file){

    ///  @brief Read the mesh from the hdf5 file
    ///  @param cmlist array containing the number of neighbours for each particle  
    ///  @param node_nbr array containing the list of neighbours for each particle  
    ///  @param input_file File name from which co-ordinate will be read
    /// 

    hid_t   file_id, dataset_id;  /* identifiers */
    herr_t  status;
    if(access(input_file.c_str(),F_OK)!=0){
        fprintf(stderr, "The configuration file does not exit\n");
        exit(1);
    }

  /* Open an existing file. */
  file_id = H5Fopen(input_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT); 

  dataset_id = H5Dopen(file_id, "cumu_list", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, 
          H5S_ALL, H5S_ALL, H5P_DEFAULT, cmlist);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "node_nbr", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, 
          H5S_ALL, H5S_ALL, H5P_DEFAULT, node_nbr);
  status = H5Dclose(dataset_id);

  status = H5Fclose(file_id);
  if(status != 0){
      fprintf(stderr, "file close failed\n");
  }
}

void io_read_config(double *Pos, 
        int N, char *file ){

    ///  @brief Read position from the file; 
    /// @note The dump should be in binary
    /// 
    FILE *fid;

    fid = fopen(file, "rb");
    if(fread(Pos, N*sizeof(double), 1, fid) != 1);
    fclose(fid);
}

void hdf5_io_read_int(int *data, string input_file, string dset_name){

    hid_t   file_id, dataset_id, space_id;  /* identifiers */
    herr_t  status;
    hsize_t dims; 
    /* Open an existing file. */
    file_id = H5Fopen(input_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT); 
    dataset_id = H5Dopen(file_id, dset_name.c_str(), H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    // Get the dataspace and check dimensions
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);
    if(status != 0){
        fprintf(stderr, "file close failed\n");
    }
    cout << "Read " << dset_name << " from " << input_file << endl;
}

// void hdf5_io_dump_int(int *stick, int N, string input_file, string dset_name){

//     hid_t   file_id, dset1, space_id;  /* identifiers */
//     herr_t  status;
//     hsize_t          dims; 

//   /* Open an existing file. */
//     dims = N;
//     file_id = H5Fcreate (input_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

//     space_id = H5Screate_simple (1, &dims, NULL);
//     dset1 = H5Dcreate2(file_id, dset_name.c_str(), H5T_NATIVE_INT, space_id, H5P_DEFAULT,
//                 H5P_DEFAULT, H5P_DEFAULT);
//     status = H5Dwrite (dset1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
//                 stick);
//     status = H5Dclose (dset1);
//     status = H5Sclose (space_id);
//     status = H5Fclose (file_id);
//   if(status != 0){
//       fprintf(stderr, "file close failed\n");
//   }
// }

// void hdf5_io_dump_bool(bool *stick, int N, 
//         string input_file, string dset_name){


//     hid_t   file_id, dset1, space_id;  /* identifiers */
//     herr_t  status;
//     hsize_t          dims; 

//   /* Open an existing file. */
//     dims = N;
//     file_id = H5Fcreate (input_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

//     space_id = H5Screate_simple (1, &dims, NULL);
//     dset1 = H5Dcreate2(file_id, dset_name.c_str(), H5T_NATIVE_HBOOL, space_id, H5P_DEFAULT,
//                 H5P_DEFAULT, H5P_DEFAULT);
//     status = H5Dwrite (dset1, H5T_NATIVE_HBOOL, H5S_ALL, H5S_ALL, H5P_DEFAULT,
//                 stick);
//     status = H5Dclose (dset1);
//     status = H5Sclose (space_id);
//     status = H5Fclose (file_id);
//   if(status != 0){
//       fprintf(stderr, "file close failed\n");
//   }
// }

void hdf5_io_read_bool(bool *stick, 
        string input_file, string dset_name){

    hid_t   file_id, dataset_id, space_id;  /* identifiers */
    herr_t  status;

    /* Open an existing file. */
    file_id = H5Fopen(input_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT); 

    dataset_id = H5Dopen(file_id, dset_name.c_str(), H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_HBOOL, 
            H5S_ALL, H5S_ALL, H5P_DEFAULT, stick);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);
    if(status != 0){
        fprintf(stderr, "file close failed\n");
    }

}

void io_dump_config_ascii(double *Pos, 
        int N, string file ){
    /// @brief dump position to the file in ascii; 
    ///
    FILE *fid;
    int i;
    fid = fopen(file.c_str(), "wb");
    for(i=0;i<N;i=i+3){
        fprintf(fid,"%g %g %g\n", Pos[i], Pos[i+1], Pos[i+3]);
    }
    fclose(fid);
}

void hdf5_io_delete(string filename){
    if (std::remove(filename.c_str()) == 0) {
        std::cout << "File deleted successfully: " << filename << std::endl;
    }
}

bool hdf5_io_has_dataset(string input_file, string dset_name) {
    ///  @brief Check if a dataset exists in the HDF5 file
    ///  @param input_file File name to check
    ///  @param dset_name Name of the dataset to check for
    ///  @return true if dataset exists, false otherwise

    if(access(input_file.c_str(), F_OK) != 0) {
        return false;  // File doesn't exist
    }

    hid_t file_id = H5Fopen(input_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        return false;  // Couldn't open file
    }

    // Check if dataset exists
    htri_t exists = H5Lexists(file_id, dset_name.c_str(), H5P_DEFAULT);
    
    H5Fclose(file_id);
    
    return exists > 0;
}