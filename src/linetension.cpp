#include "linetension.hpp"
#include "misc.hpp"
using namespace std;
extern "C" void LineTensionRead(double *, char *);
LTN::LTN(const MESH_p& mesh, string fname): mesh(mesh){
    char tmp_fname[128];
    string parafile = fname+"/para_file.in";
    sprintf(tmp_fname, "%s", parafile.c_str() );
    LineTensionRead(&lambda, tmp_fname);
    
    ofstream out_;
    out_.open( fname+"/linetensionpara.out");
    out_<< "# =========== Line tension parameters ==========" << endl
        << " N " << mesh.N << endl
        << " lambda " << lambda << endl;
    out_.close();
};

// Version 1: compute interface size

double LTN::energy_ipart(double *lijsq, int idx){
    double lt=0e0;
    int ghost=mesh.nghst;
    bool lipA=mesh.compA[idx];
    int jdx;
    for (int j = 0; j < mesh.numnbr[idx]; ++j){
        jdx = mesh.node_nbr_list[idx*ghost+j];
        if (lipA!=mesh.compA[jdx]){
            lt+=sqrt(lijsq[j]);
        }
    }
    return lt*lambda;
}

double LTN::energy_ipart(int idx){
    int ghost=mesh.nghst;
    bool lipA=mesh.compA[idx];
    int jdx;
    int num_nbr = mesh.numnbr[idx];
    double lijsq[num_nbr];
    int j;
    Vec3d rij;
    if (mesh.btype[idx] == PBC){
        for (int i =0; i < num_nbr; i++){
            j = mesh.node_nbr_list[idx*ghost+i];
            rij = diff_pbc(mesh.pos[idx], mesh.pos[j], mesh.boxlen);
            lijsq[i] = inner_product(rij, rij);
        }
    }else{
        for (int i =0; i < num_nbr; i++){
            j = mesh.node_nbr_list[idx*ghost+i];
            rij = mesh.pos[idx] - mesh.pos[j];
            lijsq[i] = inner_product(rij, rij);
        }
    }
    return energy_ipart(lijsq,idx);
}

// Version 2: count number of interfaces assigning each edge a unit length

// double LTN::energy_ipart(int idx){
//     bool lipA=mesh.compA[idx];
//     int jdx;
//     int num_nbr = mesh.numnbr[idx];
//     double lijsq[num_nbr];
//     int j;
//     Vec3d rij;
//     int ghost=mesh.nghst;
//     double lt=0e0;
//     for (int j = 0; j < mesh.numnbr[idx]; ++j){
//         jdx = mesh.node_nbr_list[idx*ghost+j];
//         if (lipA!=mesh.compA[jdx]){
//             lt+=1;
//         }
//     }
//     return lambda*lt;
// }

// double LTN::energy_total(){
//     double lt_tot=0;
//     for (int i = 0; i < mesh.N; ++i){
//         lt_tot+=energy_ipart(i);
//     }
//     return lt_tot/2;
// }

double LTN::energy_total(){
    double lt_tot=0;
    for (int i = 0; i < mesh.N; ++i){
        lt_tot+=energy_ipart(i);
    }
    return lt_tot/2;
}