#include <hdf5.h>
#include "../include/global.h"
#include "misc.h"
  /* The input for the hdf5 configuration */ 
  /* supposedly will read all the position and */ 
  /* triangles info from a single file */

void hdf5_io_read_pos(double *Pos, int *cmlist,
        int *node_nbr, int2 *bond_nbr, int *triangles,
        char input_file[]){

    hid_t   file_id, group_id,dataset_id;  /* identifiers */
    herr_t  status;

  /* Open an existing file. */
  file_id = H5Fopen(input_file, H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset_id = H5Dopen(file_id, "pos", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, 
          H5S_ALL, H5S_ALL, H5P_DEFAULT,Pos);
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);
}


void hdf5_io_read_config(double *Pos, int *cmlist,
        int *node_nbr, int2 *bond_nbr, int *triangles,
        string input_file){

    hid_t   file_id, group_id,dataset_id;  /* identifiers */
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
}

int io_dump_config(POSITION *Pos, double len, 
        int iter, int N){
    char part_file[80];
    int i;
    double x_n, y_n;
    FILE *fid;


    sprintf(part_file,"output/part_%04d.dat",iter);

    fid = fopen(part_file, "wb");

    for(i=0; i < N; i++){
        x_n = fmod((Pos[i].x + 30*len), len);
        y_n = fmod((Pos[i].y + 30*len), len);
        fprintf(fid, "%lf  %lf \n", x_n, y_n);
        fflush(fid);
    }
    fclose(fid);
    return 1;
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

}
