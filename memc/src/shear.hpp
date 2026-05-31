#ifndef SHEAR_HPP
#define SHEAR_HPP
#include <string>
#include "vector.hpp"
#include "mesh.hpp"

class SHEAR {
public:
    int initSHEAR(int N, std::string outfolder);
    void shear_positions(Vec3d *Pos, int N);
    double frame_spring_energy_ipart(Vec3d pos, Vec3d pos_t0);
    bool do_shear();
private:
    bool _do_shear;
    double constant;
    double slope;
    int shear_every;
};

#endif
