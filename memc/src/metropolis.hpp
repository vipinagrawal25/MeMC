#ifndef METROPOLIS_HPP
#define METROPOLIS_HPP
#include <string>
#include <fstream>
#include "mesh.hpp"
#include "vector.hpp"
#include "bending.hpp"
#include "stretching.hpp"
#include "sticking.hpp"
#include "activity.hpp"
#include "mdcelllist.hpp"
#include "shear.hpp"

class McP {
public :
  McP (BE &beobj, STE &steobj, STICK &stickobj,
       ACT &actobj, MDCellList &celllistobj, SHEAR &shearobj): beobj(beobj), steobj(steobj), stickobj(stickobj),
                                              actobj(actobj), celllistobj(celllistobj), shearobj(shearobj) {};
  int monte_carlo_3d(Vec3d *pos, MESH_p mesh);
  double energy_mc_3d(Vec3d *pos, MESH_p mesh,  int );
  int monte_carlo_fluid(Vec3d *, MESH_p, double);
  bool Boltzman(double DE, double activity);
  bool Glauber(double DE, double activity);
  int initMC(int, std::string);
  bool isfluid();
  bool isrestart();
  int fluidizeevery();
  int dumpskip();
  int totaliter();
  int onemciter();
  double evalEnergy(Vec3d *, MESH_p, fstream &, int);
  double getarea();
  double getvolume();
  void setEneVol(double );
  bool issemisolid();
  int* getsolidIdx();
  void mark_solid_neighbours(MESH_p mesh);
private:
    BE &beobj;
    STE &steobj;
    STICK &stickobj;
    ACT &actobj;
    MDCellList &celllistobj;
    SHEAR &shearobj;
    std::string algo;
    double dfac;
    int one_mc_iter, tot_mc_iter, dump_skip;
    double kBT;
    double delta; // increment of position
    bool is_restart;
    bool is_fluid;
    bool is_semisolid;
    int num_solid_points;
    int *solid_idx;
    int min_allowed_nbr;
    int fluidize_every;
    double fac_len_vertices;
    double totEner, totvol; 
    double EneMonitored, VolMonitored;
    double volt0;
    int acceptedmoves;
};

#endif
