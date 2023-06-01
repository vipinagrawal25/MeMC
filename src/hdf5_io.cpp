#include <hdf5.h>
#include "global.h"
#include "misc.h"
#include <unistd.h>
 /**  
 *  @brief hdf5 IO for the mesh  
 *  
 */
 
void hdf5_io_write_pos(double *Pos, int N, string input_file){

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
    dset1 = H5Dcreate2(file_id, "/pos", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT,
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

void hdf5_io_read_pos(double *Pos, string input_file){

    ///  @brief Read from the hdf5 file
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param input_file File name from which co-ordinate will be read
    ///
    hid_t   file_id,dataset_id;  /* identifiers */
    herr_t  status;

    if(access(input_file.c_str(),F_OK)!=0){
        cout << input_file.c_str() << endl;
        fprintf(stderr, "The configuration file does not exit\n");
        exit(1);
    }
  /* Open an existing file. */
  file_id = H5Fopen(input_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  dataset_id = H5Dopen(file_id, "pos", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, 
          H5S_ALL, H5S_ALL, H5P_DEFAULT,Pos);
  status = H5Dclose(dataset_id);
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

// void io_dump_config(double *Pos, 
//         int N, char *file ){
//     FILE *fid;

//     ///  @brief dump position to the file; 
//     /// @note The dump is in binary
//     /// 

//     fid = fopen(file, "wb");
//     fwrite(Pos, N*sizeof(double), 1, fid);
//     fclose(fid);
// }

// void io_dump_config_ascii(double *Pos, 
//         int N, char *file ){

//     /// @brief dump position to the file in ascii; 
//     /// 

//     FILE *fid;
//     int i;
//     fid = fopen(file, "wb");
//     for(i=0;i<N;i=i+2){
//         fprintf(fid,"%g %g\n", Pos[i], Pos[i+1]);
//     }
//     fclose(fid);
// }