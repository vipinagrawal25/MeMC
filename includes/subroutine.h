#ifndef subroutine_h
#define subroutine_h
#include "Vector.h"
#include "global.h"
#include "misc.h"
#include <iostream>
#include <random>
#include <iomanip>
#include <sstream>
using namespace std;
//*************************************************//
// metropolis.cpp
int monte_carlo_3d(Vec3d *pos, Vec3d *pos_t0, MESH_p mesh, 
                double *lij_t0, MBRANE_p mbrane,
                MC_p mcpara, AREA_p,  STICK_p ,  VOL_p , AFM_p afm, 
                ACTIVE_p activity,  SHEAR_p shear);

 double energy_mc_3d(Vec3d *pos, Vec3d *, MESH_p mesh, 
         double *lij_t0, int idx, MBRANE_p , STICK_p ,
         VOL_p , AFM_p , SHEAR_p );
int monte_carlo_surf2d(Vec2d *Pos, 
        Nbh_list *neib, LJ_p para, 
        MC_p mcpara, char *metric);

int monte_carlo_fluid(Vec3d *, MESH_p , MBRANE_p , MC_p, FLUID_p );

void init_rng(uint32_t seed_val);

//*************************************************//
//forces_lj.c
void make_nlist(Vec2d *Pos, Nbh_list *neib,
        LJ_p para, char *metric);
bool len_check(Vec2d s1, Vec2d s2, double len, 
        double new_rc, char* metric, int);

double pairlj_ipart_energy(Vec2d *Pos, int *n_list,
        int ni, int i_p, LJ_p para, char *metric);

double pairlj_total_energy(Vec2d *Pos, Nbh_list *neib,
        LJ_p para, char *metric);
//forces_surf.c
double vol_energy_change(MBRANE_p mbrane, VOL_p , double dvol);
double bending_energy_total(Vec3d *pos, MESH_p mesh, MBRANE_p para);
Vec2d bending_energy_ipart(Vec3d *pos, int *node_nbr,  
        int num_nbr, int idx, MBRANE_p para);
double bending_energy_ipart_neighbour(Vec3d *pos, 
        MESH_p mesh, int idx, MBRANE_p para);

void identify_attractive_part(Vec3d *pos, 
        bool *is_attractive, double theta_attr, int N);

double stretch_energy_total(Vec3d *pos,
         MESH_p mesh, double *lij_t0, MBRANE_p,  AREA_p );

double stretch_energy_ipart(Vec3d *pos,
         int *node_nbr, double *lij_t0, int num_nbr,
                             int idx, AREA_p para);

double stick_bottom_surface(Vec3d, Vec3d, STICK_p); 
double lj_bottom_surface(double zz, bool is_attractive, 
        double sur_pos, double eps, double sigma);
double volume_total(Vec3d *pos, 
        MESH_p mesh, MBRANE_p para);
double lj_bottom_surf_total(Vec3d *pos, Vec3d *, 
         MBRANE_p para, STICK_p st_p);
double volume_ipart(Vec3d *pos, int *node_nbr, 
        int num_nbr, int idx);
double lj_afm(Vec3d , AFM_p);
double lj_afm_total(Vec3d *pos, Vec3d *afm_force,
        MBRANE_p para, AFM_p afm);
double area_total(Vec3d *, MESH_p ,
         MBRANE_p );
double area_ipart(Vec3d *, int *, int , int);
      
double frame_spring_energy(Vec3d , Vec3d, SHEAR_p );
double spring_tot_frame_ener(Vec3d *, MBRANE_p , SHEAR_p );
int monte_carlo_shear(Vec3d *pos,  Vec3d *, MESH_p , double *, 
        MBRANE_p mbrane, ACTIVE_p activity,
        MC_p mcpara, AREA_p, SHEAR_p shear);


//init.c
//
void shear_positions(Vec3d *Pos, int N, SHEAR_p shear );
void init_system_random_pos(Vec2d *Pos,  double len, int N, char *metric, int);
double PV_change(double ,double );
/* double spring_energy(Vec3d pos, int idx, MESH_p mesh, SHEAR_p shear); */
/* double spring_tot_energy_force(Vec3d *Pos, Vec3d *spring_force, */ 
                               /* MESH_p mesh, SHEAR_p shear); */
//initialise.c
void init_eval_lij_t0(Vec3d *Pos, MESH_p mesh,
         double *lij_t0, MBRANE_p *para, SHEAR_p *shear, bool );
void init_stick_bottom(Vec3d *, MESH_p , STICK_p , 
        FLUID_p , MBRANE_p );

void init_read_config();
void init_afm_tip(AFM_p );
void init_read_parameters(MBRANE_p *mbrane_para, MC_p *mc_para, AREA_p *, FLUID_p *fld_para, 
        VOL_p *vol_para, STICK_p *stick_para, AFM_p *afm_para,  ACTIVE_p *act_para, 
        SHEAR_p *shear_para, string para_file);
 
void init_activity(ACTIVE_p, int );
int randint(int n);
void write_parameters(MBRANE_p mbrane, MC_p mc_para, AREA_p , FLUID_p fld_para, 
        VOL_p vol_p, STICK_p stick_para, AFM_p afm_para,  ACTIVE_p act_para, 
        SHEAR_p shear_para, string out_file);
 
//hdf5_io
void hdf5_io_write_pos(double *Pos, int N, string input_file);
void hdf5_io_read_pos(double *Pos, string input_file);
void hdf5_io_read_mesh(int *cmlist, int *node_nbr, string input_file);
void hdf5_io_write_mesh(int *cmlist,
        int *node_nbr, int N, int ng, string output_file);
void hdf5_io_dump_stickidx(int *stick, int N, string input_file);
void io_dump_config(double *Pos, int N, char *);
void io_read_config(double *Pos, int N, char *);
void io_dump_config_ascii(double *Pos, int N, char *);
//cubic_solve
int cubic_solve(double, double, double, double, double *);
#endif
