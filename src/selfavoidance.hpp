#ifndef SELFAVOIDANCE_HPP
#define SELFAVOIDANCE_HPP

#include <vector>
#include <array>
#include <cmath>
#include "mesh.hpp"
#include "vector.hpp"
#include "modules/celllist.hpp"
// A structure representing a particle in 3D space
// struct Particle {
//     double x, y, z;  // Position
//     std::array<double, 3> velocity; // Velocity in 3D
//     std::array<double, 3> force;    // Force in 3D

//     Particle(double x, double y, double z) : x(x), y(y), z(z), velocity({0.0, 0.0, 0.0}), force({0.0, 0.0, 0.0}) {}
// };

// A class to handle the 3D cell list algorithm for molecular dynamics
class SelfAvoid: private CellList {
public:
  SelfAvoid(string fname);
  double computeSelfRep(MESH_p , int );
  double totalRepulsiveEnergy(MESH_p);
  bool isSelfRepulsive() {return doselfrepulsion;}
  private:
  double sig, epsl;
  bool doselfrepulsion;
};
#endif // MDCELLLIST_HPP