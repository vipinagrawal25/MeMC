#include <hdf5.h>
#include "global.h"
#include "misc.h"
  /* The input for the hdf5 configuration */ 
  /* supposedly will read all the position and */ 
  /* triangles info from a single file */

void hdf5_io_read_pos(double *Pos, int *cmlist,
        int *node_nbr, int2 *bond_nbr, int *triangles,
        char input_file[]){

    hid_t   file_id,dataset_id;  /* identifiers */
    herr_t  status;

  /* Open an existing file. */
  file_id = H5Fopen(input_file, H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset_id = H5Dopen(file_id, "pos", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, 
          H5S_ALL, H5S_ALL, H5P_DEFAULT,Pos);
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);

  if(status != 0){
      fprintf(stderr, "file close failed");
  }
}


void hdf5_io_read_config(double *Pos, int *cmlist,
        int *node_nbr, int2 *bond_nbr, int *triangles,
        string input_file){

    hid_t   file_id,dataset_id;  /* identifiers */
    herr_t  status;
  /* Open an existing file. */
  if (FileExists(input_file)){
    file_id = H5Fopen(input_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT); 
  }else{
    cout << "# HDF5 file is not found to start the simulaton." << endl;
    cout << "# EXITING the code" << endl;
    exit(1);
  }
  dataset_id = H5Dopen(file_id, "pos", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, 
          H5S_ALL, H5S_ALL, H5P_DEFAULT,Pos);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "cumu_list", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, 
          H5S_ALL, H5S_ALL, H5P_DEFAULT, cmlist);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "node_nbr", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, 
          H5S_ALL, H5S_ALL, H5P_DEFAULT, node_nbr);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "bond_nbr", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, 
          H5S_ALL, H5S_ALL, H5P_DEFAULT, bond_nbr);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "triangles", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, 
          H5S_ALL, H5S_ALL, H5P_DEFAULT, triangles);
  status = H5Dclose(dataset_id);

  status = H5Fclose(file_id);
  if(status != 0){
      fprintf(stderr, "file close failed");
  }
}

void io_read_config(double *Pos, 
        int N, char *file ){
    FILE *fid;

    fid = fopen(file, "rb");
    fread(Pos, N*sizeof(double), 1, fid);
    fclose(fid);
}

void io_dump_config(double *Pos, 
        int N, char *file ){
    FILE *fid;


    fid = fopen(file, "wb");
    fwrite(Pos, N*sizeof(double), 1, fid);
    fclose(fid);
}

void io_dump_config_ascii(double *Pos, 
        int N, char *file ){
    FILE *fid;
    int i;
    fid = fopen(file, "wb");
    for(i=0;i<N;i=i+2){
        fprintf(fid,"%g %g\n", Pos[i], Pos[i+1]);
    }
    fclose(fid);
}

void hdf5_io_dump_restart_config(double *Pos, int *cmlist,
        int *node_nbr, int2 *bond_nbr, 
        int *triangles, MBRANE_para mbrane,
        string folder){

    string filename;
    string syscmds;
    int err; 
    hid_t                   file_id, dset1, space_id;
    herr_t                  status;
    hsize_t                 dims[1]; 

    syscmds="mv "+folder+"/restart.h5 "+folder+"/restart_old.h5";
    cout << syscmds << endl;
    err = system(syscmds.c_str());
    if(err !=0 ) fprintf(stderr, "fail to execute system cmds");
    filename=folder+"/restart.h5";
    /*
     * Create a new file using the default properties.
     */
    dims[0] = 3*mbrane.N;
    /* dims[1] = 3*mbrane.N; */
    file_id = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    space_id = H5Screate_simple (1, dims, NULL);
    dset1 = H5Dcreate2(file_id, "/pos", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                Pos);
    status = H5Dclose (dset1);
    status = H5Sclose (space_id);

    dims[0] = mbrane.N+1;
    space_id = H5Screate_simple (1, dims, NULL);
    dset1 = H5Dcreate2(file_id, "cumu_list", H5T_NATIVE_INT, space_id, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                cmlist);
    status = H5Dclose (dset1);
    status = H5Sclose (space_id);


    dims[0] = mbrane.num_nbr;
    space_id = H5Screate_simple (1, dims, NULL);
    dset1 = H5Dcreate2(file_id, "node_nbr", H5T_NATIVE_INT, space_id, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                node_nbr);
    status = H5Dclose (dset1);
    status = H5Sclose (space_id);


    dims[0] = 2*mbrane.num_nbr;
    space_id = H5Screate_simple (1, dims, NULL);
    dset1 = H5Dcreate2(file_id, "bond_nbr", H5T_NATIVE_INT, space_id, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                bond_nbr);
    status = H5Dclose (dset1);
    status = H5Sclose (space_id);


    dims[0] = 3*mbrane.num_triangles;
    space_id = H5Screate_simple (1, dims, NULL);
    dset1 = H5Dcreate2(file_id, "triangles", H5T_NATIVE_INT, space_id, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                triangles);
    status = H5Dclose (dset1);
    status = H5Sclose (space_id);


    status = H5Fclose (file_id);

    if(status != 0){
        fprintf(stderr, "file close failed");
    }


}
