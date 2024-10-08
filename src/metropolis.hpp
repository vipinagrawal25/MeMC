#ifndef METROPOLIS_HPP
#define METROPOLIS_HPP
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <functional>

#include "mesh.hpp"
#include "vector.hpp"
#include "bending.hpp"
#include "stretching.hpp"
#include "multicomp.hpp"
#include "electrostatics.hpp"

using namespace std;
class McP{
public : 
  McP (BE &beobj, STE &steobj, MulCom &lipidobj, ESP &chargeobj);
  McP (BE &beobj, STE &steobj);
  int monte_carlo_3d(Vec3d *pos, MESH_p mesh);
  int monte_carlo_fluid(Vec3d *, MESH_p);
  int monte_carlo_lipid(Vec3d *pos, MESH_p mesh);
  int initMC(MESH_p, std::string);
  bool isfluid();
  bool isrestart();
  int fluidizeevery();
  int dumpskip();
  int totaliter();
  int onemciter();
  double evalEnergy(MESH_p mesh);
  void write_energy(fstream &, int, const MESH_p&);
  double getarea();
  double getvolume(); 
  void setEneVol();
private:
  BE &beobj;
  STE &steobj;
  MulCom &lipidobj;
  ESP &chargeobj;
  std::string algo;
  double dfac;
  int one_mc_iter, tot_mc_iter, dump_skip;
  double kBT;
  double delta; // increment of position
  bool is_restart;
  bool is_fluid;
  int min_allowed_nbr;
  int fluidize_every;
  double fac_len_vertices;
  double totEner, totvol, bende, stretche, pre, regsole=0, electroe, vole=0;
  double EneMonitored, VolMonitored;
  double volt0;
  int acceptedmoves;
  bool sphere;
  function<double(vector<double> &, Vec3d*, MESH_p, int)> energy_mc_3d;
  double energy_mc_best(vector<double>& energy, Vec3d *pos, MESH_p mesh, int idx);
  double energy_mc_bestch(vector<double>& energy, Vec3d *pos, MESH_p mesh, int idx);
  function<bool(double, double)> Algo;
  bool Boltzman(double DE, double activity);
  bool Glauber(double DE, double activity);
};
#endif