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
int monte_carlo_3d(Vec3d *pos, MESH_p mesh, 
                double *lij_t0, MBRANE_p mbrane,
                MC_p mcpara, STICK_p ,  VOL_p , AREA_p, AFM_p afm, 
                ACTIVE_p activity,  SPRING_p spring, SPCURV_p spcurv);
double energy_mc_3d(Vec3d *pos, MESH_p mesh, 
         double *lij_t0, int idx, MBRANE_p , STICK_p ,
         VOL_p , AFM_p , SPRING_p , SPCURV_p );
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
double bending_energy_total(Vec3d *pos, MESH_p mesh, MBRANE_p para, SPCURV_p spcurv_para);
double bending_energy_ipart(Vec3d *pos, int *node_nbr,  
        int num_nbr, int idx, MBRANE_p para, SPCURV_p spcurv_para);
double bending_energy_ipart_neighbour(Vec3d *pos, 
        MESH_p mesh, int idx, MBRANE_p para, SPCURV_p spcurv_para);
void identify_attractive_part(Vec3d *pos, 
        bool *is_attractive, double theta_attr, int N);
double stretch_energy_total(Vec3d *pos,
         MESH_p mesh, double *lij_t0, MBRANE_p para);
double stretch_energy_ipart(Vec3d *pos,
         int *node_nbr, double *lij_t0, int num_nbr,
                             int idx, MBRANE_p para);
double voronoi_area(double cotJ, double cotK, 
        double jsq, double ksq, double area);
double lj_bottom_surface(double zz,
        bool is_attractive, 
        double sur_pos, double eps, double sigma);
double volume_total(Vec3d *pos, 
        MESH_p mesh, MBRANE_p para);
double lj_bottom_surf_total(Vec3d *pos,
         MBRANE_p para, STICK_p st_p);
double volume_ipart(Vec3d *pos, int *node_nbr, 
        int num_nbr, int idx);
double lj_afm(Vec3d , AFM_p);
double lj_afm_total(Vec3d *pos, Vec3d *afm_force,
        MBRANE_p para, AFM_p afm);
double area_energy_ipart(Vec3d *pos, int *node_nbr, double *area_t0, int num_nbr,
                        int idx, double coef_area_expansion);
void area_ipart(Vec3d *pos, double *area, int *node_nbr, int num_nbr, int idx);
//init.c
void init_system_random_pos(Vec2d *Pos,  double len, int N, char *metric, int);
double PV_change(double ,double );
double spring_energy(Vec3d pos, int idx, MESH_p mesh, SPRING_p spring);
double spring_tot_energy_force(Vec3d *Pos, Vec3d *spring_force, 
                               MESH_p mesh, SPRING_p spring);
void init_spcurv(SPCURV_p spcurv, Vec3d *pos, int N);
void init_area_t0(Vec3d *pos, MESH_p mesh, MBRANE_p mbrane_para, AREA_p area_para);
double area_energy_total(Vec3d *pos, MESH_p mesh, MBRANE_p para, AREA_p area_para);
//initialise.c
void init_eval_lij_t0(Vec3d *Pos, MESH_p mesh,
         double *lij_t0, MBRANE_p *para, SPRING_p *spring, bool);
void init_read_config();
void init_afm_tip(AFM_p );
bool init_read_parameters(MBRANE_p *mbrane_para, SPCURV_p *spcurv_para, MC_p *mc_para,
        FLUID_p *fld_para, VOL_p *vol_para, AREA_p *area_para, STICK_p *stick_para, 
        AFM_p *afm_para, ACTIVE_p *act_para, SPRING_p *spring_para, string para_file);
void init_activity(ACTIVE_p, int );
int randint(int n);
void write_parameters(MBRANE_p mbrane, SPCURV_p spcurv_para, MC_p mc_para, 
        FLUID_p fld_para, VOL_p vol_p, AREA_p area_p, STICK_p stick_para, AFM_p afm_para,  
        ACTIVE_p act_para, SPRING_p spring_para, string out_file);
//hdf5_io
void hdf5_io_write_pos(double *Pos, int N, string input_file);
void hdf5_io_read_pos(double *Pos, string input_file);
void hdf5_io_read_mesh(int *cmlist, int *node_nbr, string input_file);
void hdf5_io_write_mesh(int *cmlist, 
    int *node_nbr, int N, int ng, string output_file);
void io_dump_config(double *Pos, int N, char *);
void io_read_config(double *Pos, int N, char *);
void io_dump_config_ascii(double *Pos, int N, char *);
//cubic_solve
int cubic_solve(double, double, double, double, double *);
#endif