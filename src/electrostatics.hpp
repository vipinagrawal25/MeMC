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
  bool calculate(){return ical;}
  bool isexch(){return iex;}
  double getch1(){return charge1;}
  double getch2(){return charge2;}
private:
  bool ical=0;
  bool iex=0;
  void initcharges(int *compA, int N);
  // double debye_huckel(const Vec3d p1, const Vec3d p2, double q1, double q2);
  std::function<double(const Vec3d, const Vec3d, double, double)> debye_huckel;
  double dh_pbc(const Vec3d p1, const Vec3d p2, double q1, double q2);
  // double dh_nopbc(const Vec3d p1, const Vec3d p2, double q1, double q2);
  double dh_channelX(const Vec3d p1, const Vec3d p2, double q1, double q2);
  double dh_channelY(const Vec3d p1, const Vec3d p2, double q1, double q2);
  double dh_NONE(const Vec3d p1, const Vec3d p2, double q1, double q2);
  double charge1, charge2;
  std::vector<double> charges;
  double conc;  // Concentration of the electrolyte
  // Debye length
  double debyelen ;
  double kappa ;
  const double lb = 0.71; // Bjerrum length
  double boxlen = 0.0;
  bool PBCx = false;
  bool PBCy = false;
  // Detect and set periodic boundary conditions in x/y from mesh
  void determinePBC(const MESH_p& mesh);
};
#endif