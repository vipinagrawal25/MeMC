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
  void startcycle(int cycle);
  void updateparam(int anneal, string fname);
  bool exchange(){return iexch;}
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
  int one_mc_iter, tot_mc_iter, dump_skip, nexch_iter;
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
  function<double(vector<double> &, Vec3d*, MESH_p, 
          int, int,
          int, int,
          int, int)> energy_mc_exch;
  function<double(vector<double> &, Vec3d*, MESH_p, int, int, int)> energy_mc_3d;
  double energy_mc_be(vector<double>& , Vec3d *, MESH_p ,
                      int , int, 
                      int, int,
                      int, int);
  double energy_mc_best(vector<double>& , Vec3d *, MESH_p , int , int, int);
  double energy_mc_bestch(vector<double>&, Vec3d *, MESH_p , int , int, int);
  double energy_mc_bestchli(vector<double>&, Vec3d *, MESH_p , int , int, int);
  double energy_mc_bestchrep(vector<double>&, Vec3d *, MESH_p , int , int, int);
  double energy_mc_bech(vector<double>&, Vec3d *, MESH_p, 
                        int, int,
                        int, int,
                        int, int);
  double energy_mc_bechli(vector<double>&, Vec3d *, MESH_p, 
                        int, int,
                        int, int,
                        int, int);
  double energy_mc_ch(vector<double>&, Vec3d *, MESH_p, 
                    int, int, 
                    int, int,
                    int, int);
  function<bool(double, double)> Algo;
  bool Boltzman(double DE, double activity);
  bool Glauber(double DE, double activity);
  void changeparam(double dfac, double kBT, bool is_restart, int tot_mc_iter, 
      int dumpskip);
  double ini_tot_mc_iter;
  string exchtype="Global";
  function<int(int, int*, int, int, int)> get_idx2;
  int local_idx(int num_nbr, int *node_nbr_list, int cm_idx);
  int global_idx(int nframe, int N);
};
#endif