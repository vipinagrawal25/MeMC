#include "electrostatics.hpp"
#include <omp.h>
#include <string>
#include <fstream>
#include "misc.hpp"
#define tiny 1e-16
extern "C" void  ElectroRead(double*, double*, double*, char *);

void ESP::initcharges(int *compA, int N){
    for (int i = 0; i < N; ++i){
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

    if (charge1||charge2) ical=1;
    else ical=0;
    if (charge1 != charge2) iex=1;
    else iex=0;

    initcharges((int *) mesh.compA, mesh.N);
    debyelen = 0.304/sqrt(conc+tiny);
    kappa = 1/debyelen;
    
    // Store mesh properties for PBC calculations
    boxlen = mesh.boxlen;
    pbc = mesh.pbc;
    
    // Assign function pointer based on boundary conditions
    if (pbc) {
        debye_huckel = [this](const Vec3d p1, const Vec3d p2, double q1, double q2) -> double {
            return this->dh_pbc(p1, p2, q1, q2);
        };
    } else {
        debye_huckel = [this](const Vec3d p1, const Vec3d p2, double q1, double q2) -> double {
            return this->dh_nopbc(p1, p2, q1, q2);
        };
    }

    std::ofstream out_;
    out_.open(fname+"/electrostatpara.out");
    out_<< "# =========== electrostatic parameters ==========" << endl
      << " N " << mesh.N << endl
      << " charge1 = " << charge1 << endl
      << " charge2 = " << charge2 << endl
      << " electrolyte conc = " << conc << endl
      << " Debye length = " << debyelen << " nm"<< endl
      << " PBC = " << (pbc ? "true" : "false") << endl
      << " box length = " << boxlen << endl;
    out_.close();
}

// Debye-Hückel potential with periodic boundary conditions
double ESP::dh_pbc(const Vec3d p1, const Vec3d p2, double q1, double q2){
    if (q1 == 0.0 || q2 == 0.0) return 0.0;
    // Use minimum image convention for periodic boundaries
    Vec3d dr = diff_pbc(p1, p2, boxlen);
    double r = sqrt(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);
    if (r == 0) return 0.0;  // Avoid self-interaction
    return (q1 * q2 / r) * exp(-kappa * r);
}

// Debye-Hückel potential without periodic boundary conditions
double ESP::dh_nopbc(const Vec3d p1, const Vec3d p2, double q1, double q2){
    if (q1 == 0.0 || q2 == 0.0) return 0.0;
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
    if (charge1 == 0.0) return 0.0;
    Vec3d Pos1 = Pos[idx];
    for (int j = 0; j < N; ++j) {
        total_potential += debye_huckel(Pos1, Pos[j], charge1, charges[j]);
    }
    return total_potential*lb;
}

double ESP::debye_huckel_total(Vec3d *Pos, int N){
    double total_potential = 0.0;
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j){
            total_potential += debye_huckel(Pos[i], Pos[j], charges[i], charges[j]);
        }
    }
    return total_potential;
}

void ESP::exchange(int idx1, int idx2){
    swap(charges[idx1], charges[idx2]);
}