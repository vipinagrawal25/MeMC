#include "electrostatics.hpp"
#include <omp.h>
#include <string>
#include <fstream>
#include "misc.hpp"

extern "C" void  ElectroRead(double*, double*, double*, char *);

void ESP::initcharges(int *compA, int N){
    // print(compA,N);
    // std::cout << N << "\n";
    for (int i = 0; i < N; ++i){
        // std::cout << i << "\n";
        if (compA[i]) {charges.push_back(charge2);}
        else {charges.push_back(charge1);}
    }
}

ESP::ESP(const MESH_p& mesh, std::string fname){
    char tmp_fname[128];
    string parafile, outfile;

    parafile = fname+"/para_file.in";
    sprintf(tmp_fname, "%s", parafile.c_str());
    ElectroRead(&charge1, &charge2, &conc, tmp_fname);

    if (charge1*charge2) ical=1;

    // print(mesh.compA,mesh.N);
    initcharges((int *) mesh.compA, mesh.N);
    
    debyelen = 0.304/sqrt(conc);
    kappa = 10/debyelen;
    lb = lb/10;

    std::ofstream out_;
    out_.open(fname+"/electrostatpara.out");
    out_<< "# =========== bending parameters ==========" << endl
      << " N " << mesh.N << endl
      << " charge1 = " << charge1 << endl
      << " charge2 = " << charge2 << endl
      << " electrolyte conc = " << conc << endl
      << " Debye length = " << debyelen << " nm"<< endl;
    out_.close();
}

double ESP::debye_huckel(const Vec3d p1, const Vec3d p2, double q1,
    double q2){
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;
    double r = sqrt(dx * dx + dy * dy + dz * dz);
    if (r == 0) return 0.0;  // Avoid self-interaction
    return (q1 * q2 / r) * exp(-kappa * r);
}

double ESP::debye_huckel_ipart(Vec3d *Pos, int idx, int N){
    double total_potential = 0.0;
    double charge1 = charges[idx];
    double charge2;
    Vec3d Pos1 = Pos[idx];
    // omp_set_num_threads(1);
    // #pragma omp parallel for reduction(+:total_potential) schedule(dynamic)
    for (int j = 0; j < N; ++j) {
        total_potential += debye_huckel(Pos1, Pos[j], charge1, charges[j]);
    }
    return total_potential*lb;
}

double ESP::debye_huckel_total(Vec3d *Pos, int N){
    double total_potential = 0.0;
    // omp_set_num_threads(1);
    // #pragma omp parallel for reduction(+:total_potential) schedule(dynamic)
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j){
            total_potential += debye_huckel(Pos[i], Pos[j], charges[i], charges[j]);
        }
    }
    return total_potential;
}

void ESP::exchange(int idx1, int idx2){
   double temp = charges[idx1];
   charges[idx1] = charges[idx2];
   charges[idx2] = temp;
}
