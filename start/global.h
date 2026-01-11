// this is a global file
#ifndef GLOBAL_H
#define GLOBAL_H
//
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include <random>
#include <fstream>
#include <mpi.h>
using namespace std;
#define pi 3.14159265358979
#define R_del 0.05
#define char_len 128
typedef struct{
    double x, y;
}Vec2d;

typedef struct{
    /// @brief MCpara Structure
    /// @param dfac weight of the random increment
    /// @param one_mc_iter number of moves such that each particle is selected
    /// atleast once 
    /// @param tot_mc_iter Number of one_mc_iter 
    /// @param dump_skip the config file is dumped every this many iters 
    /// @param kBT Boltzmann  constant times temperature
    /// @param delta approximated value of lattice spacing
    double dfac;
    int one_mc_iter;
    int tot_mc_iter;
    int dump_skip;
    double kBT;
    double delta; // increment of position
    bool is_restart;
}MC_p;


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
#endif
