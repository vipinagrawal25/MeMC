// this is a global file
#ifndef GLOBAL_H
#define GLOBAL_H
//
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
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
 }MCpara;

typedef struct{
    string act;
    double minA, maxA;
    double *activity;
}ActivePara;

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
    /// @param num_triangles;  //number of triangles 
    /// @param num_nbr; //  neighbours of all particles
    double coef_bend;  //coefficient bending
    double YY;  //coefficient stretching
    double coef_vol_expansion;   //coefficient of volume expansion
    double pressure;
    double radius;  // radius of ball
    double sigma, epsilon, theta; // sigma and epsilon for the bottom attractive wall
    double sp_curv; // spontaneous curvature of the membrane.
    double pos_bot_wall;  // position of the bottom attractive wall
    double av_bond_len; // average length of the bond
    double *tot_energy;
    double *volume; // these are updated after each monte carlo 
    int N;   // number of particles in mesh
    int num_triangles;  //number of triangles 
    int num_nbr; //  neighbours of all particles
    bool istick;
}MBRANE_para;
//

//
typedef struct{
    /// @brief Mesh Structure
    /// @param size_node_nbr Number of neighbours 2*N - 4 
    /// @param cmlist; Cumulative sum of neighbours 
    /// @param node_nbr_list; list of neighbours of a node 
    int size_node_nbr;
    int *cmlist;
    int *node_nbr_list;
    int nPole, sPole;
}MESH;
//

typedef struct{
    /// @brief Afm structure
    /// @param tip_pos_z;  position of tip in z
    /// @param tip_rad;  radius of the tip
    /// @param sigma, epsilon; 
    /// sigma epsilon of the LJ potential used to model AFM
    double tip_pos_z; // position of tip in z
    double tip_rad; // radius of the tip
    double sigma, epsilon;
    int icompute;
}AFM_para;
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
}LJpara;
//
typedef struct{
    int cnt;
    int list[200];
}Nbh_list;
//
typedef struct{
    int icompute;
    double constant;
    double nPole_eq_z;
    double sPole_eq_z;
}SPRING_para;
#endif