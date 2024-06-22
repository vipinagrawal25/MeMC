#ifndef BENDING_HPP
#define BENDING_HPP
#include <string>
#include "Mesh.hpp"
#include "Vector.hpp"
// #include "global.h"

class BE {
public : 
  double bending_energy_ipart_neighbour(Vec3d *pos, MESH_p mesh, int idx);
  double bending_energy_ipart(Vec3d *pos, int *node_nbr, int num_nbr,
                              int idx);
  double bending_energy_total(Vec3d *pos, MESH_p mesh);
  int initBE(int, std::string);
private:
    double coef_bend;
    string which_spcurv;
    double minC, maxC, theta;
    double spcurv;
};

#endif
