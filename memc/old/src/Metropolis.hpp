#ifndef METROPOLIS_H
#define METROPOLIS_H
#include "global.h"

extern "C" void  MC_listread(char *, double *, double *, bool *,
           int *, int *, char *);
class MC_P{
public : 
int monte_carlo_3d(Vec3d *pos, MESH_p mesh, 
                double *lij_t0, MBRANE_p mbrane,
                STICK_p ,  VOL_p, AREA_p, AFM_p afm, 
                ACTIVE_p activity,  SPRING_p spring, SPCURV_p spcurv);
double energy_mc_3d(Vec3d *pos, MESH_p mesh, 
         double *lij_t0, int idx, MBRANE_p , STICK_p ,
         VOL_p , AFM_p , SPRING_p , SPCURV_p );
int monte_carlo_fluid(Vec3d *, MESH_p , MBRANE_p , MC_p, FLUID_p );
int initMC(int, string);
private:
    string algo;
    double dfac;
    int one_mc_iter, tot_mc_iter; dump_skip;
    double kBT;
    double delta; // increment of position
};

#endif
