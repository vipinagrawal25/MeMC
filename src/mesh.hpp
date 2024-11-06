#ifndef MESH_HPP
#define MESH_HPP
#include <string>
#include "vector.hpp"
#include "hdf5_io.hpp"
#include <vector>
#include "misc.hpp"
#include <cmath>

extern "C" void MeshRead(int *, int *, int *, double *, char*, char *);

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
    std::string distribution;
    double radius, ini_vol;
    MESH_p(std::string outfolder){
        char tmp_fname[128], tmp_dist[128];;
        std::string para_file = outfolder+"/para_file.in";
        sprintf(tmp_fname, "%s", para_file.c_str());

        MeshRead(&bdry_type, &nghst, &ncomp, &compfrac, tmp_dist, tmp_fname);
        distribution=tmp_dist;
        N = (int)hdf5_io_get_Np(outfolder+"/input.h5", "pos")/3;

        pos = new Vec3d[N];
        numnbr = new int[N];
        node_nbr_list = new int[N*nghst];
        compA = new int[N];

        hdf5_io_read_double( (double *)pos,  outfolder+"/input.h5", "pos" );
        hdf5_io_read_mesh((int *) numnbr, (int *)node_nbr_list, 
            outfolder+"/input.h5");

        if (isPlaner()){
            sphere=false;
            edge = get_nstart(N, 1);
            boxlen=get_box_dim().first*(1+1/sqrt(N));
        }
        else{
            sphere=true;
            edge = -1;
            boxlen=0;
            bdry_type=2;   // Sphere always has a pbc.
            radius=calculateRadius();
            ini_vol = 4e0/3e0*M_PI*radius*radius*radius;
            cout << radius << endl;
        }

        if (ncomp==1) compfrac=0;
        
        if (distribution=="Random" || distribution=="random"){
            fillPoints(compA, compfrac, N);    
        }else if(distribution=="Janus" || distribution=="janus"){
            identify_attractive_part(compA, pos, 2*M_PI*compfrac, N);
        }

    }

    void free(){
        delete[] pos;
        delete[] node_nbr_list;
        delete[] numnbr;
        delete[] compA;
    }

    pair<double, double> get_box_dim(){
        // Vec3d *Pos = mesh.pos;
        double xmin=1e7,xmax=-1e7,ymin=1e7,ymax=-1e7;
        for (int i = 0; i < N; ++i){
            if (pos[i].x<xmin)  xmin=pos[i].x;
            if (pos[i].x>xmax)  xmax=pos[i].x;
            if (pos[i].y<ymin)  ymin=pos[i].y;
            if (pos[i].y>ymax)  ymax=pos[i].y;
        }
        double xlen = xmax-xmin;
        double ylen = ymax-ymin;
    return {xlen, ylen};
    }

    bool isPlaner(){
        for(int i=0;i<N;i++){if(pos[i].z!=0){return false;}}
    return true;
    }

    double calculateRadius() {
        double sum_distances = 0.0;
    
        for (int i=0; i<N; i++) {
            auto point=pos[i];
            double distance = std::sqrt(point.x * point.x + point.y * point.y + point.z * point.z);
            sum_distances += distance;
        }
    
        // Calculate average distance (approximate radius)
        double radius = sum_distances / N;
    return radius;
    }

};
#endif