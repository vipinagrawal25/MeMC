#ifndef subroutine_h
#define subroutine_h
#include "Position.h"
#include "../src/misc.h"
#include <iostream>
//*************************************************//
// main.c
int monte_carlo_3d(POSITION *pos, MESH mesh, 
                double *lij_t0, bool *, MBRANE_para mbrane, 
                MCpara mcpara);
//*************************************************//          
//forces_lj.c
void make_nlist();
void make_nlist_pf();
bool len_check_ss(POSITION s1, POSITION s2, double len, double new_rc);
bool len_check_pf(POSITION s1, POSITION s2, double len, double new_rc, char* metric);
double cal_length();
double pairlj_ipart_energy(POSITION *Pos, int *n_list,
        int ni, int i_p, LJpara para, char *metric);
double pairlj_total_energy();
double pairlj_ipart_energy_pf(POSITION *Pos, 
        int i_p, LJpara para, char *metric);
double pairlj_total_energy_pf(POSITION *Pos, LJpara para, char *metric);
/* double pairlj(); */

//forces_surf.c
double bending_energy_total(POSITION *pos, MESH mesh, MBRANE_para para);
double bending_energy_ipart(POSITION *pos, int *node_nbr, int2 *bond_nbr, 
        int num_nbr, int idx, MBRANE_para para,
        string method="curvature");
double bending_energy_ipart_neighbour(POSITION *pos, 
        MESH mesh, int idx, MBRANE_para para);
double stretch_energy_total(POSITION *pos, 
        MESH mesh, double *lij_t0,
         MBRANE_para para);
void identify_obtuse(POSITION *pos, int *triangles, 
       double *obtuse,  int N);
void identify_attractive_part(POSITION *pos, 
        bool *is_attractive, int N);
double stretch_energy_ipart(POSITION *pos, 
        int *node_nbr, double *lij_t0,
        int num_nbr, int idx, MBRANE_para para);
double lj_bottom_surface(double zz,
        bool is_attractive, 
        double sur_pos, double eps, double sigma);

double lj_bottom_surf_total(POSITION *pos, 
         bool *is_attractive, MBRANE_para para);

void volume_area_enclosed_membrane(POSITION *pos, 
    int *triangles, int num_triangles,
    double *avolume, double *aarea);

double volume_ipart(POSITION *pos, 
        int *node_nbr, int2* bond_nbr,
        int num_nbr, int idx, MBRANE_para para);

double lj_afm(POSITION , AFM_para);
double lj_afm_total(POSITION *pos, POSITION *afm_force,
        MBRANE_para para, AFM_para afm);

//initialise.c
void initialize_system();
void initialize_eval_lij_t0(POSITION *Pos, MESH mesh, 
        double *lij_t0, MBRANE_para *para);
void initialize_read_config();
void initialize_afm_tip(AFM_para );
void initialize_read_parameters( MBRANE_para *mbrane, 
        AFM_para *afm, MCpara *mcpara, string para_file);
int randint(int n);

//visit_io.c
void visit_vtk_io(double *points, 
        int *triangles, 
        int Np, string filename, 
        string dataname="miketesting");
void visit_vtk_io_point_data(bool *data, 
        int Np, string filename, 
        string dataname);
void visit_vtk_io_cell_data(double *data, 
        int Np, string filename, 
        string dataname);

void visit_vtk_io_afm_tip(double *data, 
        int Np, string filename);
//hdf5_io
void hdf5_io_read_pos(double *Pos, int *cmlist,
        int *node_nbr, int2 *bond_nbr, int *triangles,
        char input_file[]);
void hdf5_io_read_config(double *Pos, int *cmlist,
        int *node_nbr, int2 *bond_nbr, int *triangles,
        string input_file);
int io_dump_config(POSITION *Pos, double len, 
        int iter, int N);
void hdf5_io_dump_restart_config(double *Pos, int *cmlist,
        int *node_nbr, int2 *bond_nbr, 
        int *triangles, MBRANE_para mbrane,
        string folder);
//misc
bool FileExists(const std::string &s);
int cubic_solve(double, double, double, double, double *);

#endif
