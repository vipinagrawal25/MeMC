#ifndef subroutine_h
#define subroutine_h
#include "global.h"
#include <iostream>
#include <random>
#include <iomanip>
#include <sstream>
using namespace std;

int monte_carlo_surf2d(Vec2d *Pos, 
        Nbh_list *neib, LJ_p para, 
        MC_p mcpara, char *metric);
void init_rng(uint32_t seed_val);
void init_system_random_pos(Vec2d *Pos,  double len, 
        int N, char *metric, int bdry_condt );
 
//forces_lj.c
void make_nlist(Vec2d *Pos, Nbh_list *neib,
        LJ_p para, char *metric);
bool len_check(Vec2d s1, Vec2d s2, double len, 
        double new_rc, char* metric, int);

double pairlj_ipart_energy(Vec2d *Pos, int *n_list,
        int ni, int i_p, LJ_p para, char *metric);

double pairlj_total_energy(Vec2d *Pos, Nbh_list *neib,
        LJ_p para, char *metric);

//initialise.c
void hdf5_io_write_pos(double *Pos, int N, string input_file);
void hdf5_io_read_pos(double *Pos, string input_file);
void io_dump_config(double *Pos, int N, char *);
void io_read_config(double *Pos, int N, char *);
void io_dump_config_ascii(double *Pos, int N, char *);
#endif
