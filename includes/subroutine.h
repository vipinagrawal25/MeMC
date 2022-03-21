#ifndef subroutine_h
#define subroutine_h
#include "Vector.h"
#include "misc.h"
#include <iostream>

//*************************************************//
// Metropolis.c
int monte_carlo_3d(Vec3d *pos, MESH mesh, 
                double *lij_t0, bool *is_attractive, 
                MBRANE_para mbrane, 
                MCpara mcpara, AFM_para afm);
int monte_carlo_surf2d(Vec2d *Pos, 
        Nbh_list *neib, LJpara para, 
        MCpara mcpara, char *metric);
void init_rng(uint32_t seed_val);

//*************************************************//          
//forces_lj.c
void make_nlist(Vec2d *Pos, Nbh_list *neib,
        LJpara para, char *metric);
bool len_check(Vec2d s1, Vec2d s2, double len, double new_rc, char* metric);
double pairlj_ipart_energy(Vec2d *Pos, int *n_list,
        int ni, int i_p, LJpara para, char *metric);
double pairlj_total_energy(Vec2d *Pos, Nbh_list *neib,
        LJpara para, char *metric);

//*************************************************//          
//forces_surf.c
double bending_energy_total(Vec3d *pos, MESH mesh, MBRANE_para para);
double bending_energy_ipart(Vec3d *pos, int *node_nbr,  
        int num_nbr, int idx, MBRANE_para para);
double bending_energy_ipart_neighbour(Vec3d *pos, 
        MESH mesh, int idx, MBRANE_para para);
double stretch_energy_total(Vec3d *pos, 
        MESH mesh, double *lij_t0,
         MBRANE_para para);
void identify_attractive_part(Vec3d *pos, 
        bool *is_attractive, double theta_attr, int N);
double stretch_energy_ipart(Vec3d *pos, 
        int *node_nbr, double *lij_t0,
        int num_nbr, int idx, MBRANE_para para);
double lj_bottom_surface(double zz,
        bool is_attractive, 
        double sur_pos, double eps, double sigma);
double volume_total(Vec3d *pos, 
        MESH mesh, MBRANE_para para);
double lj_bottom_surf_total(Vec3d *pos, 
         bool *is_attractive, MBRANE_para para);
double volume_ipart(Vec3d *pos, 
        int *node_nbr, 
        int num_nbr, int idx, MBRANE_para para);

double lj_afm(Vec3d , AFM_para);
double lj_afm_total(Vec3d *pos, Vec3d *afm_force,
        MBRANE_para para, AFM_para afm);

//init.c
void init_system_random_pos(Vec2d *Pos,  double len, int N, char *metric);
void init_eval_lij_t0(Vec3d *Pos, MESH mesh, 
        double *lij_t0, MBRANE_para *para);
void init_read_config();
void init_afm_tip(AFM_para );

void init_read_parameters( MBRANE_para *mbrane, 
        AFM_para *afm, MCpara *mcpara, char *para_file);
int randint(int n);

//hdf5_io
void hdf5_io_write_pos(double *Pos, int N, char *input_file);
void hdf5_io_read_pos(double *Pos, char *input_file);
void hdf5_io_read_mesh(int *cmlist,
        int *node_nbr, char *input_file);
void io_dump_config(double *Pos, int N, char *);
void io_read_config(double *Pos, int N, char *);
void io_dump_config_ascii(double *Pos, int N, char *);

//cubic_solve
int cubic_solve(double, double, double, double, double *);

#endif
