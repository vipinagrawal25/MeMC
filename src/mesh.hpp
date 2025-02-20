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

extern "C" void MeshRead(int *, int *, int *, double *, char*, char *);
using namespace std;
struct MESH_p{
    /// @brief Mesh Structure
    /// @param numnbr; number of neighbours
    /// @param node_nbr_list; list of neighbours of a node
    int N, bdry_type;
    int nghst;
    bool sphere;
    int *numnbr;
    int *node_nbr_list;
    double boxlen;
    int lastbdry; // storing corner index specially for periodic case.
    double av_bond_len;
    Vec3d *pos;
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
    MESH_p(string outfolder){
        char tmp_fname[128], tmp_dist[128];;
        string para_file = outfolder+"/para_file.in";
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
            auto allbonds = make_bond_list(node_nbr_list);
            // Currently, the code saves the last index of the boundary, assuming that the boundary indices are at the beginning.
            lastbdry = get_bdry(allbonds);
            boxlen=get_box_dim().first*(1+1/sqrt(N));
            pbc = determine_pbc();
            if (pbc){
                if(bdry_type == 0 || bdry_type == 1) {
                    cerr << "Error: Boundary type of fixed frame (0) or channel (1) cannot be used with periodic mesh. Run this code with the right boundary condition" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }   
        else{
            sphere=true;
            lastbdry = -1;
            boxlen=0;
            bdry_type=2;   // Sphere always has a free boundary condition.
            radius=calculateRadius();
            cout << "radius = " << radius << endl;
            ini_vol = 4e0/3e0*M_PI*radius*radius*radius;
        }

        if (ncomp==1) compfrac=0;
        if (compfrac==0) ncomp=1;
        if (distribution=="Random" || distribution=="random"){
            fillPoints(compA, compfrac, N);
        }else if(distribution=="Janus" || distribution=="janus"){
            zattr=2*(compfrac-0.5) * radius;
            identify_attractive_part(compA, pos, zattr, N);
        }else if (distribution=="Point" || distribution=="point"){
            compA[0]=1;
        }

        for(i = 0; i < N; i++){
            num_nbr = numnbr[i];
            cm_idx = nghst * i;
            for(k = cm_idx; k < cm_idx + num_nbr; k++){
                j = node_nbr_list[k];
                dr = diff_pbc(pos[j], pos[i], boxlen);
                sum_lij += sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
                npairs++;
            }
        }
        av_bond_len = sum_lij/npairs;
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
    
    double calculateRadius(){
        double sum_distances = 0.0;
    
        for (int i=0; i<N; i++) {
            auto point=pos[i];
            double distance = sqrt(point.x * point.x + point.y * point.y + point.z * point.z);
            sum_distances += distance;
        }
    
        // Calculate average distance (approximate radius)
        double radius = sum_distances / N;
    return radius;
    }
    
    set<pair<int,int>> make_bond_list(int *node_nbr){
        set<pair<int,int>> edge_set;
        // Determine the number of vertices from the length of neighbor_indices and nghst.    
        for (int i = 0; i < N; ++i) {
            int current_index = i * nghst;
            // Loop over the valid neighbors for vertex i
            for (int j = 0; j < numnbr[i]; ++j) {
                int neighbor = node_nbr[current_index + j];
                // Create the edge (i, neighbor) in sorted order
                pair<int,int> edge = make_pair(min(i, neighbor), max(i, neighbor));
                edge_set.insert(edge);
            }
        }
        return edge_set;
    }

    int get_bdry(const set<pair<int,int>>& edge_set){
        for (int i = 0; i < N; ++i){
            int* nbrs = node_nbr_list + i * nghst;
            vector<int> nbr_list(nbrs, nbrs + numnbr[i]);
            if (is_boundary_vertex(nbr_list, edge_set) == false) return i-1;
        }
        return -1;
    }

    bool is_boundary_vertex(const vector<int>& nbrs, const set<pair<int,int>>& edge_set){
        size_t n = nbrs.size();
        if (n < 3) {
            return true;
        }
        for (size_t i = 0; i < n; ++i) {
            int n1 = nbrs[i];
            int n2 = nbrs[(i + 1) % n];
            // Create an edge with sorted order
            pair<int,int> edge = make_pair(min(n1, n2), max(n1, n2));
            if (edge_set.find(edge) == edge_set.end()) {
                return true;
            }
        }
        return false;
    }

    bool determine_pbc(){
        // Check if any node has a neighbor that wraps around the boundary
        for (int i = 0; i < N; ++i) {
            int current_index = i * nghst;
            for (int j = 0; j < numnbr[i]; ++j) {
                int neighbor = node_nbr_list[current_index + j];
                if (abs(pos[i].x - pos[neighbor].x) > boxlen / 2 ||
                    abs(pos[i].y - pos[neighbor].y) > boxlen / 2 ||
                    abs(pos[i].z - pos[neighbor].z) > boxlen / 2) {
                    return true;
                }
            }
        }
        return false;
    }

};

#endif
