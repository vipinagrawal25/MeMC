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
double bending_energy_ipart(Vec3d *pos, int *node_nbr, int2 *bond_nbr, 
        int num_nbr, int idx, MBRANE_para para,
        string method="curvature");
double bending_energy_ipart_neighbour(Vec3d *pos, 
        MESH mesh, int idx, MBRANE_para para);
double stretch_energy_total(Vec3d *pos, 
        MESH mesh, double *lij_t0,
         MBRANE_para para);
void identify_obtuse(Vec3d *pos, int *triangles, 
       double *obtuse,  int N);
void identify_attractive_part(Vec3d *pos, 
        bool *is_attractive, int N);
double stretch_energy_ipart(Vec3d *pos, 
        int *node_nbr, double *lij_t0,
        int num_nbr, int idx, MBRANE_para para);
double lj_bottom_surface(double zz,
        bool is_attractive, 
        double sur_pos, double eps, double sigma);

double lj_bottom_surf_total(Vec3d *pos, 
         bool *is_attractive, MBRANE_para para);

void volume_area_enclosed_membrane(Vec3d *pos, 
    int *triangles, int num_triangles,
    double *avolume, double *aarea);

double volume_ipart(Vec3d *pos, 
        int *node_nbr, int2* bond_nbr,
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
        AFM_para *afm, MCpara *mcpara, string para_file);
int randint(int n);

//vtk_io.c
void vtk_io(double *points, 
        int *triangles, 
        int Np, string filename, 
        string dataname="miketesting");
void vtk_io_point_data(bool *data, 
        int Np, string filename, 
        string dataname);
void vtk_io_cell_data(double *data, 
        int Np, string filename, 
        string dataname);

//hdf5_io
void hdf5_io_read_pos(double *Pos, int *cmlist,
        int *node_nbr, int2 *bond_nbr, int *triangles,
        char input_file[]);
void hdf5_io_read_config(double *Pos, int *cmlist,
        int *node_nbr, int2 *bond_nbr, int *triangles,
        string input_file);
void io_dump_config(double *Pos, int N, char *);
void io_read_config(double *Pos, int N, char *);
void io_dump_config_ascii(double *Pos, int N, char *);
void hdf5_io_dump_restart_config(double *Pos, int *cmlist,
        int *node_nbr, int2 *bond_nbr, 
        int *triangles, MBRANE_para mbrane,
        string folder);
//misc
bool FileExists(const std::string &s);
int cubic_solve(double, double, double, double, double *);

#endif
