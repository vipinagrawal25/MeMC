#ifndef MESH_HPP
#define MESH_HPP

#include <string>
#include "vector.hpp"
#include "hdf5_io.hpp"
#include <vector>
#include "misc.hpp"
#include <cmath>
#include <iostream>
#include <set>
#include <algorithm>
#include <utility>  // for pair, make_pair

extern "C" void MeshRead(int *, double *, char*, char *);

using namespace std;
struct MESH_p{
    /// @brief Mesh Structure
    /// @param numnbr; number of neighbours
    /// @param node_nbr_list; list of neighbours of a node

    int N, bdry_type;
    const int nghst=12;
    bool sphere;
    int *numnbr;
    int *node_nbr_list;
    double boxlen;
    int lastbdry; // storing corner index specially for periodic case.
    double av_bond_len;
    Vec3d *pos;
    int *cells;
    int ncomp;
    double compfrac;
    int *compA;
    string distribution;
    double radius, ini_vol, zattr;
    double sum_lij = 0.0;
    int npairs = 0;
    Vec3d dr;
    int num_nbr, cm_idx, i, j, k;
    bool pbc;
    /// Declare all the functions here.
    bool isPlaner();
    // set<pair<int, int>> make_bond_list(int *node_nbr_list);
    // int get_bdry(const set<pair<int, int>> &allbonds);
    pair<double, double> get_box_dim();
    bool determine_pbc();
    double calculateRadius();    
    // Constructor
    MESH_p(string outfolder);
    // Cleanup method
    void free();
};

#endif