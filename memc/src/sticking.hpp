#ifndef STICKING_HPP
#define STICKING_HPP
#include <string>
#include <vector>
#include "mesh.hpp"
#include "vector.hpp"

// #include "global.h"

class STICK {
public : 
  double stick_energy_total(Vec3d *pos, int);
 double stick_energy_ipart(Vec3d pos);
 int initSTICK(int , std::string );
private:
    double eps1, eps2;  // strength of stick potential
    double sigma, pos_bot_wall;
};

#endif
