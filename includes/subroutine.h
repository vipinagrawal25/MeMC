#ifndef subroutine_h
#define subroutine_h
#include "Vector.h"
#include "misc.h"
#include <iostream>
#include <random>
#include <iomanip>
#include <sstream>
using namespace std;
//*************************************************//
// metropolis.cpp
int monte_carlo_3d(Vec3d *pos, MESH mesh, 
                double *lij_t0, bool *is_attractive, 
                MBRANE_para mbrane,
                MCpara mcpara, AFM_para afm, ActivePara activity,  SPRING_para spring);
bool Metropolis(double DE, double, MCpara mcpara);
// double rand_inc_theta(double th0, double dfac);
 double energy_mc_3d(Vec3d *pos, MESH mesh, 
         double *lij_t0, bool *is_attractive, int idx, bool *is_be_pos,
         MBRANE_para mbrane,MCpara mcpara, AFM_para afm,
         SPRING_para spring);
int monte_carlo_surf2d(Vec2d *Pos, 
        Nbh_list *neib, LJpara para, 
        MCpara mcpara, char *metric);

int monte_carlo_fluid(Vec3d *, MESH , MBRANE_para , MCpara , AFM_para ,
                ActivePara , SPRING_para );

void init_rng(uint32_t seed_val);

//*************************************************//
//forces_lj.c
void make_nlist(Vec2d *Pos, Nbh_list *neib,
        LJpara para, char *metric);
bool len_check(Vec2d s1, Vec2d s2, double len, 
        double new_rc, char* metric, int);

double pairlj_ipart_energy(Vec2d *Pos, int *n_list,
        int ni, int i_p, LJpara para, char *metric);

double pairlj_total_energy(Vec2d *Pos, Nbh_list *neib,
        LJpara para, char *metric);
//forces_surf.c
double vol_energy_change(MBRANE_para mbrane,double dvol);
double bending_energy_total(Vec3d *pos, MESH mesh, MBRANE_para para);
double bending_energy_ipart(Vec3d *pos, int *node_nbr,  
        int num_nbr, int idx, MBRANE_para para);
double bending_energy_ipart_neighbour(Vec3d *pos, 
        MESH mesh, int idx, MBRANE_para para);

void identify_attractive_part(Vec3d *pos, 
        bool *is_attractive, double theta_attr, int N);

 double stretch_energy_total(Vec3d *pos,
         MESH mesh, double *lij_t0, MBRANE_para para);

 double stretch_energy_ipart(Vec3d *pos,
         int *node_nbr, double *lij_t0, int num_nbr,
                             int idx, MBRANE_para para);

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
void init_system_random_pos(Vec2d *Pos,  double len, int N, char *metric, int);
double PV_change(MBRANE_para mbrane,double dvol);
double spring_energy(Vec3d pos, int idx, MESH mesh, SPRING_para spring);
double spring_tot_energy_force(Vec3d *Pos, Vec3d *spring_force, 
                               MESH mesh, SPRING_para spring);
//initialise.c
 void init_eval_lij_t0(Vec3d *Pos, MESH mesh,
         double *lij_t0, MBRANE_para *para, SPRING_para *spring);
void init_read_config();
void init_afm_tip(AFM_para );
void init_read_parameters( MBRANE_para *mbrane, 
        AFM_para *afm, MCpara *mcpara, ActivePara *, SPRING_para *spring, 
        string para_file);
void init_activity(ActivePara, int );
int randint(int n);
void write_param(string fname, MBRANE_para mbrane, MCpara mcpara, 
                SPRING_para spring);
//hdf5_io
void hdf5_io_write_pos(double *Pos, int N, string input_file);
void hdf5_io_read_pos(double *Pos, string input_file);
void hdf5_io_read_mesh(int *cmlist,
        int *node_nbr, string input_file);
void io_dump_config(double *Pos, int N, char *);
void io_read_config(double *Pos, int N, char *);
void io_dump_config_ascii(double *Pos, int N, char *);
//cubic_solve
int cubic_solve(double, double, double, double, double *);
#endif
