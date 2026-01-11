#ifndef SELFAVOIDANCE_HPP
#define SELFAVOIDANCE_HPP

#include <vector>
#include <array>
#include <cmath>
#include "mesh.hpp"
#include "vector.hpp"
#include "modules/sparse_celllist.hpp"
// # include "celllist.hpp" -- if you do not want sparse representation

class SelfAvoid: public CellList {
public:
  SelfAvoid(MESH_p, string fname);
  double computeSelfRep(MESH_p , int );
  double totalRepulsiveEnergy(MESH_p);
  bool isSelfRepulsive() {return doselfrepulsion;}
private:
  double LJ(Vec3d p1, Vec3d p2);
  double sig, epsl;
  bool doselfrepulsion;
};
#endif // MDCELLLIST_HPP