#ifndef MDCELLLIST_HPP
#define MDCELLLIST_HPP

#include <vector>
#include <array>
#include <cmath>
// #include <iostream>
#include "mesh.hpp"
#include "vector.hpp"

// A structure representing a particle in 3D space
// struct Particle {
//     double x, y, z;  // Position
//     std::array<double, 3> velocity; // Velocity in 3D
//     std::array<double, 3> force;    // Force in 3D

//     Particle(double x, double y, double z) : x(x), y(y), z(z), velocity({0.0, 0.0, 0.0}), force({0.0, 0.0, 0.0}) {}
// };

// A class to handle the 3D cell list algorithm for molecular dynamics
class MDCellList {
public:
  MDCellList( string fname);
  // int initCelllist(string fname);
  void buildCellList(Vec3d *, int N);
  double computeSelfRep(Vec3d *, MESH_p , int );
  void printParticlesInCells(Vec3d *);
  double totalRepulsiveEnergy(Vec3d *, MESH_p);
  bool isSelfRepulsive();
  private:
  double boxSize_;    // Size of the simulation box
  double cutoff_;     // Cutoff radius for interactions
  int numCells_;      // Number of cells per dimension
  double cellSize_;   // Size of each cell
  double minL;
  double sig, epsl;
  bool doselfrepulsion;
  std::vector<std::vector<int>> cellList_;  // Cell list storing particle indices
  // std::vector<Particle> particles_;         // List of particles
  // Helper function to get the cell index based on particle position
  int getCellIndex(double x, double y, double z, double);
  // Get neighboring cells' indices (including itself)
  std::vector<int> getNeighboringCells(int cellIdx);
};

#endif // MDCELLLIST_HPP

