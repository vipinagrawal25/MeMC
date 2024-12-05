#include "linetension.hpp"
#include "misc.hpp"
using namespace std;
extern "C" void LineTensionRead(double *, char *);
LTN::LTN(const MESH_p& mesh, string fname): mesh(mesh){
    char tmp_fname[128];
    string parafile = fname+"/para_file.in";
    sprintf(tmp_fname, "%s", parafile.c_str() );
    LineTensionRead(&lambda, tmp_fname);

    cout << "lambda = " << lambda << endl;

    ofstream out_;
    out_.open( fname+"/linetensionpara.out");
    out_<< "# =========== Line tension parameters ==========" << endl
        << " N " << mesh.N << endl
        << " lambda " << lambda << endl;
    out_.close();
    // exit(1);
};

double LTN::energy_ipart(double *lijsq, int idx){
    double lt=0e0;
    int ghost=mesh.nghst;
    bool lipA=mesh.compA[idx];
    for (int j = 0; j < mesh.numnbr[idx]; ++j){
        if (lipA!=mesh.compA[idx*ghost+j]){
            lt+=lijsq[j];            
        }
    }
    return lt*lambda;
}

double LTN::energy_ipart(int idx){
    double lijsq[mesh.numnbr[idx]];
    return line_tension_ipart(lijsq,idx);
}

double LTN::energy_total(){
    double lt_tot=0;
    for (int i = 0; i < mesh.N; ++i){
        lt_tot+=line_tension_ipart(i);
    }
    return lt_tot/2;
}