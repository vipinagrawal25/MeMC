#include "shear.hpp"
#include <cmath>
#include <fstream>
#include <cstdio>

#define pi 3.14159265358979

extern "C" void Shear_listread(bool *, int *, double *, double *, char *);

int SHEAR::initSHEAR(int N, std::string outfolder) {
    char tmp_fname[128];
    std::string parafile = outfolder + "/para_file.in";
    sprintf(tmp_fname, "%s", parafile.c_str());
    Shear_listread(&_do_shear, &shear_every, &slope, &constant, tmp_fname);
    std::ofstream out_(outfolder + "/shearpara.out");
    out_ << "# =========== Shear Parameters ===========" << std::endl
         << " do_shear = " << _do_shear  << std::endl
         << " slope = "    << slope      << std::endl
         << " constant = " << constant   << std::endl;
    out_.close();
    return 0;
}

void SHEAR::shear_positions(Vec3d *Pos, int N) {
    for (int i = 0; i < N; i++)
        Pos[i].x += slope * (Pos[i].y - pi);
}

double SHEAR::frame_spring_energy_ipart(Vec3d pos, Vec3d pos_t0) {
    double xmin = pos_t0.x + slope * (pos_t0.y - pi);
    return 0.5 * constant * pow(pos.x - xmin, 2);
}

bool SHEAR::do_shear() { return _do_shear; }
