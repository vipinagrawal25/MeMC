#ifndef ELECTROSTAT_HPP
#define ELECTROSTAT_HPP
#include <string>
#include "mesh.hpp"
#include "vector.hpp"
#include "misc.hpp"
#include <vector>
// #include "global.h"
class ESP{
public:
  ESP(const MESH_p& mesh, std::string fname);
  double debye_huckel_ipart(Vec3d *Pos, int idx, int N);
  double debye_huckel_total(Vec3d *Pos, int N);
  void exchange(int idx1, int idx2);
  bool ical=0;
private:
  void initcharges(int *compA, int N);
  double debye_huckel(const Vec3d p1, const Vec3d p2, double q1, double q2);
  double charge1, charge2;
  std::vector<double> charges;
  double conc;  // Concentration of the electrolyte
  // Debye length
  double debyelen ;
  double kappa ;
  double lb = 0.7; // Bjerrum length
};
#endif