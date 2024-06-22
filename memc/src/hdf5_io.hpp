#ifndef HDF5_IO_HPP
#define HDF5_IO_HPP
#include <hdf5.h>
#include <unistd.h>
#include "Vector.hpp"

void hdf5_io_write_double(double *Pos, int N, string input_file, string);
void hdf5_io_read_double(double *Pos, string input_file, string);
void hdf5_io_read_mesh(int *cmlist, int *node_nbr, string input_file);
void hdf5_io_write_mesh(int *cmlist,
        int *node_nbr, int N, int ng, string output_file);
void io_dump_config(double *Pos, int N, char *);
void io_read_config(double *Pos, int N, char *);
void io_dump_config_ascii(double *Pos, int N, char *);

#endif
