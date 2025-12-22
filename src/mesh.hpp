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

extern "C" void MeshRead(double *, char*, double *, char *);
enum BoundaryType { PBC, FBC, NONE };

using namespace std;
struct MESH_p{
    /// @brief Mesh Structure
    /// @param numnbr; number of neighbours
    /// @param node_nbr_list; list of neighbours of a node
    int N = 0;
    const int nghst=12;
    bool sphere = false;
    int *numnbr = nullptr;
    int *node_nbr_list = nullptr;
    double av_bond_len = 0.0;
    Vec3d *pos = nullptr;
    int *cells = nullptr;
    double compfrac = 0.0;
    int *compA = nullptr;
    string distribution;
    double radius = 0.0, ini_vol = 0.0, zattr = 0.0;
    double sum_lij = 0.0;
    int npairs = 0;
    Vec3d dr;
    int num_nbr = 0, cm_idx = 0, i = 0, j = 0, k = 0;
    bool periodic_x = false;
    bool periodic_y = false;
    int ncomp;
    bool isPlaner();
    pair<double, double> get_box_dim() const;
    double calculateRadius();
    MESH_p(string outfolder);
    BoundaryType* btype=nullptr;
    void free();
    pair<bool, bool> determinePBC();
    double boxlen ;
};

#endif