#include "global.h"
#include "subroutine.h"
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <random>
#include <sstream>
#include <sys/types.h>
#include <unistd.h>
std::mt19937 rng;
//
void init_rng(u_int32_t seed_val){
    rng.seed(seed_val);
}

int get_nstart(int N, int bdrytype){
    static int nf1;
    int nf2;
    nf1 = (int) sqrt((double)N);
    switch (bdrytype) {
        case 0:
            nf2 = 2 * nf1;
            break;
        case 1:
            nf2 = 4 * nf1; 
            break;
        default:
            nf2 = 0;
    }
    return nf2;;
}


bool Metropolis(double DE,  MC_p mcpara) {
  /// @brief Metropolis algorithm
  /// @param DE change in energy
  /// @param kbt boltzmann constant times temperature
  /// @return True if DE< 0 or the random number generated is less than
  /// exp(-DE/kbt)
  /// @details see
  /// https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm
  bool yes;
  double rand;
  std::uniform_real_distribution<> rand_real(0,1);
  yes = (DE <= 0.e0);
  if (!yes) {
    rand = rand_real(rng);
    yes = rand < exp(-DE / mcpara.kBT);
  }
  return yes;
}

double rand_inc_theta(double th0, double dfac) {
  /// @brief increment the polar angle randomly
  double dth;
  double tmp_th0;

  std::uniform_real_distribution<> rand_real(0,1);
  tmp_th0 = 10;
  while (tmp_th0 > pi || tmp_th0 < 0) {
    dth = (pi / dfac) * rand_real(rng);
    tmp_th0 = th0 + dth;
  }
  return dth;
}

//
int monte_carlo_surf2d(Vec2d *Pos, Nbh_list *neib, LJ_p para, MC_p mcpara,
                       char *metric) {
  /// @brief Monte-Carlo routine for the points on the surface of sphere or flat
  /// plane
  ///  @param Pos array containing co-ordinates of all the particles
  ///  @param neib Neighbour list information for the particle
  ///  @param para parameters related to LJ potential
  /// @param mcpara Monte-Carlo related parameters
  ///  @param metric Topology of the surface "cart" for flat plane "sph" for
  /// sphere
  /// @return number of accepted moves
  int i, move;
  double x_o, y_o, x_n, y_n;
  double de, Eini, Efin;
  double dxinc, dyinc;
  bool is_sph, is_cart;
  int n_ghost;

  n_ghost = get_nstart(para.N, para.bdry_condt);
  std::uniform_int_distribution<uint32_t> rand_int(n_ghost, para.N - 1);
  std::uniform_real_distribution<> rand_real(-1, 1);
  is_sph = false;
  is_cart = false;
  if (strcmp(metric, "sph") == 0) {
    is_sph = true;
  }
  if (strcmp(metric, "cart") == 0) {
    is_cart = true;
  }
  move = 0;
  for (i = 0; i < mcpara.one_mc_iter; i++) {
    int idx = rand_int(rng);
    Eini = pairlj_ipart_energy(Pos, neib[idx].list, neib[idx].cnt, idx, para,
                               metric);
    x_o = Pos[idx].x;
    y_o = Pos[idx].y;

    if (is_cart) {
      dxinc = (para.sigma / mcpara.dfac) * rand_real(rng);
      x_n = fmod((x_o + dxinc + 30 * para.len), para.len);
      Pos[idx].x = x_n;

      dyinc = (para.sigma / mcpara.dfac) * rand_real(rng);
      y_n = fmod((y_o + dyinc + 30 * para.len), para.len);
      Pos[idx].y = y_n;
    }
    if (is_sph) {
      dxinc = rand_inc_theta(Pos[idx].x, mcpara.dfac);
      dyinc = (para.sigma / mcpara.dfac) * (rand_real(rng));
      x_n = x_o + para.sigma * dxinc;
      y_n = fmod((y_o + dyinc + 30 * 2 * pi), 2 * pi);
      Pos[idx].x = x_n;
      Pos[idx].y = y_n;
    }

    Efin = pairlj_ipart_energy(Pos, neib[idx].list, neib[idx].cnt, idx, para,
                               metric);
    de = (Efin - Eini);
    if (Metropolis(de, mcpara)) {
      move = move + 1;
    } else {
      Pos[idx].x = x_o;
      Pos[idx].y = y_o;
    }
  }
  return move;
}
