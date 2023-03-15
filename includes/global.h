// this is a global file
#ifndef GLOBAL_H
#define GLOBAL_H
//
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <cassert>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include <cstdlib>
#include <random>
#include "Vector.h"
#include <fstream>
#include <mpi.h>
using namespace std;
#define pi 3.14159265358979
#define R_del 0.05
#define char_len 64
//not included the celid and particle id which i shall do in the cell linked list part

///  @brief Discriptions of the structure classes 
typedef struct{
    /// @brief MCpara Structure
    /// @param dfac weight of the random increment
    /// @param one_mc_iter number of moves such that each particle is selected
    /// atleast once 
    /// @param tot_mc_iter Number of one_mc_iter 
    /// @param dump_skip the config file is dumped every this many iters 
    /// @param kBT Boltzmann  constant times temperature
    /// @param delta approximated value of lattice spacing
    string algo;
    double dfac;
    int one_mc_iter;
    int tot_mc_iter;
    int dump_skip;
    double kBT;
    double delta; // increment of position
    bool is_restart;
 }MC_p;

typedef struct{
    string act;
    double minA, maxA;
    double *activity;
}ACTIVE_p;

//
typedef struct{
    /// @brief Membrane Structure
    //
    /// @param coef_bend;  coefficient bending
    /// @param coef_str;  coefficient stretching
    /// @param coef_vol_expansion;   coefficient of volume expansion
    /// @param radius;   radius of ball
    /// @param sigma, epsilon, theta;  sigma and epsilon for the bottom attractive wall
    /// @param sp_curv;  spontaneous curvature of the membrane.
    /// @param pos_bot_wall;  position of the bottom attractive wall
    /// @param av_bond_len;  average length of the bond
    /// @param *tot_energy; total energy of the system
    /// @param *volume;  these are updated after each monte carlo 
    /// @param N;   // number of particles in mesh
    ///
    double coef_bend;  //coefficient bending
    double radius;  // radius of ball
    double sp_curv; // spontaneous curvature of the membrane.
    double av_bond_len; // average length of the bond
    double *tot_energy;
    double *volume; // these are updated after each monte carlo 
    double *area; // these are updated after each monte carlo 
    int N;   // number of particles in mesh
    double len;
    int bdry_type;
}MBRANE_p;
//
//
typedef struct{
    bool is_tether;
    double YY;
    double Ka;
}AREA_p;

typedef struct{
    bool is_fluid;
    int min_allowed_nbr;
    int fluidize_every;
    double fac_len_vertices;
    // fac_len_vertices time the average length
}FLUID_p;

typedef struct{
    bool do_volume;
    bool is_pressurized;
    double coef_vol_expansion;   //coefficient of volume expansion
    double pressure;
}VOL_p;

typedef struct{
    bool do_stick;
    double pos_bot_wall;  // position of the bottom attractive wall
    double sigma, epsilon, theta; // sigma and epsilon for the bottom attractive wall
    bool *is_attractive;
}STICK_p;


//
typedef struct{
    /// @brief Mesh Structure
    /// @param numnbr; number of neighbours
    /// @param node_nbr_list; list of neighbours of a node
    int nghst;
    int *numnbr;
    int *node_nbr_list;
    int nPole, sPole;
}MESH_p;
//

typedef struct{
    /// @brief Afm structure
    /// @param tip_pos_z;  position of tip in z
    /// @param tip_rad;  radius of the tip
    /// @param sigma, epsilon; 
    /// sigma epsilon of the LJ potential used to model AFM
    bool do_afm;
    double tip_pos_z; // position of tip in z
    double tip_rad; // radius of the tip
    double sigma, epsilon;
    int icompute;
}AFM_p;
//

typedef struct{
  /// @brief LJ parameters for initial montecarlo in start.cpp
    /// @param N; number of particles 
    /// @param sigma, epsilon; sigma epsilon of the LJ potential used to model AFM
    /// @param len length of the box
    /// @param r_cut cutoff length of the potential taken as 4 sigma
    int N;
    double sigma;
    double epsilon;
    double len;
    double r_cut;
    int bdry_condt;
}LJ_p;

//
typedef struct{
    int cnt;
    int list[200];
}Nbh_list;
//
typedef struct{
    bool do_spring;
    int icompute;
    double constant;
    double nPole_eq_z;
    double sPole_eq_z;
}SPRING_p;

#endif
