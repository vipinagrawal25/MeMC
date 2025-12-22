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
#include "selfavoidance.hpp"
#include "linetension.hpp"

using namespace std;
class McP{
public : 
  McP (BE &beobj, STE &steobj, MulCom &lipidobj, ESP &chargeobj, 
    SelfAvoid &repulsiveobj, LTN &lineobj);
  int monte_carlo_3d(Vec3d *pos, MESH_p mesh);
  int monte_carlo_fluid(Vec3d *pos, MESH_p mesh);
  int monte_carlo_lipid(Vec3d *pos, MESH_p mesh);
  int initMC(MESH_p mesh, std::string fname);
  bool isfluid();
  bool isrestart();
  int fluidizeevery();
  int dumpskip();
  int totaliter();
  int onemciter();
  double evalEnergy(MESH_p mesh);
  void write_energy(fstream &fileptr, int itr, const MESH_p &mesh);
  double getarea();
  double getvolume(); 
  void setEneVol();
  void startcycle(int cycle);
  void updateparam(int anneal, string fname);
  bool exchange(){return iexch;}
  int ncycles(){return nanneal_cycle;}
  void wHeader(const MESH_p &mesh, std::fstream &fid);
  std::string initial_l0;
private:
  BE &beobj;
  STE &steobj;
  MulCom &lipidobj;
  ESP &chargeobj;
  SelfAvoid &repulsiveobj;
  LTN &lineobj;
  std::string algo;
  double dfac;
  double initial_dfac;
  int one_mc_iter, tot_mc_iter, dump_skip, nexch_iter, nanneal_cycle;
  double kBT;
  double initial_kBT;
  double delta; // increment of position
  bool is_restart, iexch;
  bool is_fluid;
  int min_allowed_nbr;
  int fluidize_every;
  double fac_len_vertices;
  double totEner, totvol, bende, stretche, pre, regsole=0, electroe, vole=0;
  double linee=0;
  double selfe=0;
  double EneMonitored, VolMonitored, AreaMonitored;
  double volt0;
  int acceptedmoves;
  bool sphere;
  std::function<double(std::vector<double>&, Vec3d*, MESH_p, int, int, int, int, int, int,
  BoundaryType, BoundaryType)> energy_mc_exch;
  std::function<double(std::vector<double>&, Vec3d*, MESH_p, int, int, int, BoundaryType)> energy_mc_3d;
  double energy_mc_be(std::vector<double>& energy, Vec3d *pos, MESH_p mesh, int idx1, int idx2, int cm_idx1, int cm_idx2, int num_nbr1, int num_nbr2, BoundaryType btype1, BoundaryType btype2);
  double energy_mc_best(std::vector<double>& energy, Vec3d *pos, MESH_p mesh, int idx, int cm_idx, int num_nbr, BoundaryType btype);
  double energy_mc_bestch(std::vector<double>& energy, Vec3d *pos, MESH_p mesh, int idx, int cm_idx, int num_nbr, BoundaryType btype);
  double energy_mc_bestchli(std::vector<double>& energy, Vec3d *pos, MESH_p mesh, int idx, int cm_idx, int num_nbr, BoundaryType btype);
  double energy_mc_bestchrep(std::vector<double>& energy, Vec3d *pos, MESH_p mesh, int idx, int cm_idx, int num_nbr, BoundaryType btype);
  double energy_mc_bech(std::vector<double>& energy, Vec3d *pos, MESH_p mesh, int idx1, int idx2, int cm_idx1, int cm_idx2, int num_nbr1, int num_nbr2, BoundaryType btype1, BoundaryType btype2);
  double energy_mc_bechli(std::vector<double>& energy, Vec3d *pos, MESH_p mesh, int idx1, int idx2, int cm_idx1, int cm_idx2, int num_nbr1, int num_nbr2, BoundaryType btype1, BoundaryType btype2);
  double energy_mc_ch(std::vector<double>& energy, Vec3d *pos, MESH_p mesh, int idx1, int idx2, int cm_idx1, int cm_idx2, int num_nbr1, int num_nbr2);
  std::function<bool(double, double)> Algo;
  bool Boltzman(double DE, double activity);
  bool Glauber(double DE, double activity);
  void changeparam(double dfac, double kBT, bool is_restart, int tot_mc_iter, int dumpskip);
  double ini_tot_mc_iter;
  std::string exchtype = "Global";
  std::function<int(int, int*, int, int, int)> get_idx2;
  int local_idx(int num_nbr, int *node_nbr_list, int cm_idx);
  int global_idx(int nframe, int N);
};
#endif