/* #include <hdf5/serial/hdf5.h> */
#include <hdf5.h>
#include "global.h"
#include "misc.h"
#include <unistd.h>
 /**  
 *  @brief hdf5 IO for the mesh  
 *  
 */

 
void hdf5_io_write_double(double *Pos, int N, 
        string input_file, string dset_name){

	 ///  @brief hdf5 io for  the position 
	 ///  @param Pos array containing co-ordinates of all the particles
	 ///  @param N number of points 
	 ///  @param input_file File name to dump all the co-ordinate
	 /// 


    hid_t   file_id, dset1, space_id;  /* identifiers */
    herr_t  status;
    hsize_t          dims; 

  /* Open an existing file. */
    dims = N;
    file_id = H5Fcreate (input_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    space_id = H5Screate_simple (1, &dims, NULL);
    dset1 = H5Dcreate2(file_id, dset_name.c_str(), H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                Pos);
    status = H5Dclose (dset1);
    status = H5Sclose (space_id);
    status = H5Fclose (file_id);
  if(status != 0){
      fprintf(stderr, "file close failed\n");
  }
}

void hdf5_io_read_double(double *Pos, string input_file, 
        string dset_name){

    ///  @brief Read from the hdf5 file
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param input_file File name from which co-ordinate will be read
    /// 


    hid_t   file_id,dataset_id;  /* identifiers */
    herr_t  status;

    if(access(input_file.c_str(),F_OK)!=0){
        fprintf(stderr, "The configuration file does not exit\n");
        exit(1);
    }
  /* Open an existing file. */
  file_id = H5Fopen(input_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  dataset_id = H5Dopen(file_id, dset_name.c_str(), H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, 
          H5S_ALL, H5S_ALL, H5P_DEFAULT,Pos);
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);

  if(status != 0){
      fprintf(stderr, "file close failed\n");
  }
}

void hdf5_io_write_mesh(int *cmlist,
        int *node_nbr, int N, int ng, string output_file){

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

    /* Open an existing file. */

    dims = N;
    space_id = H5Screate_simple (1, &dims, NULL);
    dset1 = H5Dcreate2(file_id, "/cumu_list", H5T_NATIVE_INT, space_id, H5P_DEFAULT,
            H5P_DEFAULT, H5P_DEFAULT);
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




void hdf5_io_read_mesh(int *cmlist,
        int *node_nbr,  string input_file){

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

void hdf5_io_read_int(int *stick, string input_file, string dset_name){


    hid_t   file_id, dataset_id, space_id;  /* identifiers */
    herr_t  status;
    hsize_t          dims; 

    /* Open an existing file. */
    file_id = H5Fopen(input_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT); 

    dataset_id = H5Dopen(file_id, dset_name.c_str(), H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_INT, 
            H5S_ALL, H5S_ALL, H5P_DEFAULT, stick);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);
    if(status != 0){
        fprintf(stderr, "file close failed\n");
    }

}

void hdf5_io_dump_int(int *stick, int N, string input_file, string dset_name){


    hid_t   file_id, dset1, space_id;  /* identifiers */
    herr_t  status;
    hsize_t          dims; 

  /* Open an existing file. */
    dims = N;
    file_id = H5Fcreate (input_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    space_id = H5Screate_simple (1, &dims, NULL);
    dset1 = H5Dcreate2(file_id, dset_name.c_str(), H5T_NATIVE_INT, space_id, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                stick);
    status = H5Dclose (dset1);
    status = H5Sclose (space_id);
    status = H5Fclose (file_id);
  if(status != 0){
      fprintf(stderr, "file close failed\n");
  }
}

void hdf5_io_dump_bool(bool *stick, int N, 
        string input_file, string dset_name){


    hid_t   file_id, dset1, space_id;  /* identifiers */
    herr_t  status;
    hsize_t          dims; 

  /* Open an existing file. */
    dims = N;
    file_id = H5Fcreate (input_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    space_id = H5Screate_simple (1, &dims, NULL);
    dset1 = H5Dcreate2(file_id, dset_name.c_str(), H5T_NATIVE_HBOOL, space_id, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset1, H5T_NATIVE_HBOOL, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                stick);
    status = H5Dclose (dset1);
    status = H5Sclose (space_id);
    status = H5Fclose (file_id);
  if(status != 0){
      fprintf(stderr, "file close failed\n");
  }
}

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
