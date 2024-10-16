#ifndef MESH_HPP
#define MESH_HPP
#include <string>
#include "vector.hpp"
#include "hdf5_io.hpp"
#include <vector>
#include "misc.hpp"

extern "C" void MeshRead(int *, int *, double *, int *, double *, char *);

struct MESH_p{
    /// @brief Mesh Structure
    /// @param numnbr; number of neighbours
    /// @param node_nbr_list; list of neighbours of a node
    int N, bdry_type;
    int nghst;
    bool sphere;
    // std::string topology;
    int *numnbr;
    int *node_nbr_list;
    double boxlen;
    int edge; // storing corner index specially for periodic case.
    double av_bond_len;
    Vec3d *pos;
    int ncomp;
    double compfrac;
    int *compA;
    double radius;
    double ini_vol;
    MESH_p(std::string outfolder){
        char tmp_fname[128];
        std::string para_file = outfolder+"/para_file.in";
        sprintf(tmp_fname, "%s", para_file.c_str());

        MeshRead(&bdry_type, &nghst, &radius, &ncomp, &compfrac, tmp_fname);
        
        N = (int)hdf5_io_get_Np(outfolder+"/input.h5", "pos")/3;

        pos = new Vec3d[N];
        numnbr = new int[N];
        node_nbr_list = new int[N*nghst];
        compA = new int[N];

        if (ncomp==1) compfrac=1;
        fillPoints(compA, compfrac, N);

        ini_vol = 4e0/3e0*M_PI*radius*radius*radius;
    }

    void free(){
        delete[] pos;
        delete[] node_nbr_list;
        delete[] numnbr;
        delete[] compA;
    }
};
#endif