#ifndef STICKING_HPP
#define STICKING_HPP
#include <string>
#include <vector>
#include "mesh.hpp"
#include "vector.hpp"

// #include "global.h"

class STICK {
public : 
   STICK(int , std::string );
  double stick_energy_total(Vec3d *pos, int);
  double stick_energy_ipart(Vec3d pos, int);
  void identify_attractive_part(Vec3d *);
private:
    double eps1, eps2;  // strength of stick potential
    double sigma, pos_bot_wall;
  std::vector<bool> isattractive;
};

#endif
