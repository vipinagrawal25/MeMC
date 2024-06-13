#include "global.h"
#include "subroutine.h"
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <random>
#include <sstream>
#include <unistd.h>

//
int del_nbr(int *nbrs, int numnbr, int idx) {
  // delet int idx between i1 and i2 in the nbrs list
  int new_numnbr, delete_here;
  bool logic;

  logic = false;
  delete_here = 0;

  // for(int i=0; i<numnbr+3; i++)printf("%d \n", nbrs[i]);
  // printf("\n\n");

  while (!logic) {
    logic = (nbrs[delete_here] == idx);
    ++delete_here;
  }

  memcpy(nbrs + delete_here - 1, &nbrs[delete_here],
         sizeof(int) * (numnbr - delete_here + 1));

  // for(int i=0; i<numnbr+3; i++)printf("%d \n", nbrs[i]);
  // printf("\n\n");

  return numnbr - 1;
}

int add_nbr(int *nbrs, int numnbr, int idx, int i1, int i2) {
  // add int idx between i1 and i2 in the nbrs list
  int insert_here;
  bool logic;

  logic = false;
  insert_here = 0;

  // for(int i=0; i<numnbr+3; i++)printf("%d \n", nbrs[i]);
  //     printf("\n\n");

  while (!logic) {
    logic = (nbrs[insert_here] == i1) || (nbrs[insert_here] == i2);
    ++insert_here;
  }

  logic = (nbrs[insert_here] == i1) || (nbrs[insert_here] == i2);
  if (logic) {
    memcpy(nbrs + insert_here, &nbrs[insert_here - 1],
           sizeof(int) * (numnbr - insert_here + 1));
    nbrs[insert_here] = idx;
  } else {
    insert_here = 0;
    memcpy(nbrs + insert_here, &nbrs[insert_here - 1],
           sizeof(int) * (numnbr - insert_here + 1));
    nbrs[insert_here] = idx;
  }

  // for(int i=0; i<numnbr+3; i++)printf("%d \n", nbrs[i]);
  //     printf("\n\n");

  return numnbr + 1;
}

bool Metropolis(double DE, double activity, MC_p mcpara) {
  /// @brief Metropolis algorithm
  /// @param DE change in energy
  /// @param kbt boltzmann constant times temperature
  /// @return True if DE< 0 or the random number generated is less than
  /// exp(-DE/kbt)
  /// @details see
  /// https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm
  bool yes;
  double rand;
  DE += activity;
  yes = (DE <= 0.e0);
  if (!yes) {
    rand = RandomGenerator::generateUniform(0.0,1.0);
    yes = rand < exp(-DE / mcpara.kBT);
  }
  return yes;
}
bool Glauber(double DE, double activity, MC_p mcpara) {
  /// @brief Glauber algorithm
  /// @param DE change in energy
  /// @param kbt boltzmann constant times temperature
  bool yes;
  double rand;
  DE += activity;
  rand = RandomGenerator::generateUniform(0.0,1.0);
  yes = rand < 1 / (1 + exp(DE / mcpara.kBT));
  return yes;
}

double rand_inc_theta(double th0, double dfac) {
  /// @brief increment the polar angle randomly
  double dth;
  double tmp_th0;

  tmp_th0 = 10;
  while (tmp_th0 > pi || tmp_th0 < 0) {
    dth = (pi / dfac) * RandomGenerator::generateUniform(0.0,1.0);
    tmp_th0 = th0 + dth;
  }
  return dth;
}

double energy_mc_3d(Vec3d *pos, MESH_p mesh, double *lij_t0, 
                    int idx, MBRANE_p mbrane, STICK_p st_p, 
                    VOL_p vol_p, AFM_p afm, SPRING_p spring,
                    SPCURV_p spcurv) {
  /// @brief Estimate the contribution from all the energies when a particle is
  /// moved randomly
  ///  @param Pos array containing co-ordinates of all the particles
  ///  @param mesh mesh related parameters -- connections and neighbours
  /// information
  /// @param lij_t0 initial distance between points of membrane
  /// @param is_attractive true if the zz sees the  bottom wall
  ///  @param idx index of ith particle which is moved;
  ///  @param mbrane  Membrane related parameters;
  /// @param mcpara Monte-Carlo related parameters
  /// @param AFM afm related parameter
  /// @return Change in Energy when idx particle is moved

  double E_b, E_s, E_stick, E_afm, E_spr;
  int cm_idx, num_nbr;

  E_b = 0.0;
  E_s = 0.0;
  E_stick = 0.0;
  E_afm = 0.0;
  E_spr = 0.0;

  cm_idx = mesh.nghst * idx;
  num_nbr = mesh.numnbr[idx];

  E_b = bending_energy_ipart(pos, (int *)(mesh.node_nbr_list + cm_idx), num_nbr,
                             idx, mbrane, spcurv);

  E_b += bending_energy_ipart_neighbour(pos, mesh, idx, mbrane, spcurv);

  E_s = stretch_energy_ipart(pos, (int *)(mesh.node_nbr_list + cm_idx),
                             (lij_t0 + cm_idx), num_nbr, idx, mbrane);
  if(st_p.do_stick)
  E_stick = lj_bottom_surface(pos[idx].z, st_p.is_attractive[idx],
      st_p.pos_bot_wall, st_p.epsilon, st_p.sigma); 

    if(afm.do_afm) E_afm = lj_afm(pos[idx], afm);

    if(spring.do_spring) E_spr = spring_energy(pos[idx], idx, mesh, spring);
  return E_b + E_s + E_stick + E_afm + E_spr;
}
//
int monte_carlo_3d(Vec3d *pos, MESH_p mesh, double *lij_t0, 
                   MBRANE_p mbrane, MC_p mcpara, STICK_p st_p,
                   VOL_p vol_p, AREA_p area_p, AFM_p afm,
                   ACTIVE_p activity, SPRING_p spring,
                   SPCURV_p spcurv) {
  /// @brief Monte-Carlo routine for the membrane
  ///  @param Pos array containing co-ordinates of all the particles
  ///  @param mesh mesh related parameters -- connections and neighbours
  /// information
  /// @param lij_t0 initial distance between points of membrane
  /// @param is_attractive true if the zz sees the  bottom wall
  ///  @param mbrane  Membrane related parameters;
  /// @param mcpara Monte-Carlo related parameters
  /// @param AFM afm related parameter
  /// @return number of accepted moves
  int i, move;
  int num_nbr, cm_idx;
  double x_o, y_o, z_o, x_n, y_n, z_n;
  double de, Eini, Efin;
  double dxinc, dyinc, dzinc;
  double vol_i, vol_f;
  double dvol, de_vol, ini_vol, de_pressure;
  double KAPPA;
  bool yes;
  int nframe;
  //
  nframe = get_nstart(mbrane.N, mbrane.bdry_type);

  ini_vol = (4. / 3.) * pi * pow(mbrane.radius, 3);
  KAPPA = vol_p.coef_vol_expansion;
  move = 0;
  for (i = 0; i < mcpara.one_mc_iter; i++) {
    int idx = RandomGenerator::intUniform(nframe, mbrane.N -1);
    cm_idx = idx*mesh.nghst;
    num_nbr = mesh.numnbr[idx];
    Eini = energy_mc_3d(pos, mesh, lij_t0, idx, mbrane, st_p, vol_p,
                        afm, spring, spcurv);
    if(vol_p.do_volume || vol_p.is_pressurized) vol_i = volume_ipart(pos,
            (int *) (mesh.node_nbr_list + cm_idx), num_nbr, idx);
    //
    x_o = pos[idx].x;
    y_o = pos[idx].y;
    z_o = pos[idx].z;
    //
    dxinc = (mcpara.delta / mcpara.dfac) * (RandomGenerator::generateUniform(-1.0,1.0));
    dyinc = (mcpara.delta / mcpara.dfac) * (RandomGenerator::generateUniform(-1.0,1.0));
    dzinc = (mcpara.delta / mcpara.dfac) * (RandomGenerator::generateUniform(-1.0,1.0));
    //
    x_n = x_o + dxinc;
    y_n = y_o + dyinc;
    z_n = z_o + dzinc;
    //
    pos[idx].x = x_n;
    pos[idx].y = y_n;
    pos[idx].z = z_n;
    //
    Efin = energy_mc_3d(pos, mesh, lij_t0, idx, mbrane, st_p, vol_p,
                        afm, spring, spcurv);
    de = (Efin - Eini);
    if(vol_p.do_volume){
      vol_f = volume_ipart(pos,
            (int *) (mesh.node_nbr_list + cm_idx), num_nbr, idx);
      dvol = 0.5*(vol_f - vol_i);
      de_vol = vol_energy_change(mbrane, vol_p, dvol);
      // cout << de << "\t";
      de = (Efin - Eini) + de_vol;
        // cout << de << endl;
    }
    if(vol_p.is_pressurized){
       vol_f = volume_ipart(pos,
            (int *) (mesh.node_nbr_list + cm_idx), num_nbr, idx);
        dvol=0.5*(vol_f - vol_i);
        de_pressure = PV_change(vol_p.pressure, dvol);
        // cout << de_pressure << endl;
        de = (Efin - Eini) + de_pressure;
    }
    // yes = Metropolis(de, 0.0, mcpara);
    if (mcpara.algo == "mpolis") {
      yes = Metropolis(de, activity.activity[idx], mcpara);
    } else if (mcpara.algo == "glauber") {
      yes = Glauber(de, activity.activity[idx], mcpara);
    }
    //
    if(yes) {
      move = move + 1;
      mbrane.tot_energy[0] += de;
      if(vol_p.do_volume || vol_p.is_pressurized) mbrane.volume[0] += dvol;
    } else {
      pos[idx].x = x_o;
      pos[idx].y = y_o;
      pos[idx].z = z_o;
    }
  }
  return move;
}
//

int monte_carlo_fluid(Vec3d *pos, MESH_p mesh, MBRANE_p mbrane, MC_p mcpara, FLUID_p fl_para) {

  /// @brief Monte-Carlo routine for the membrane
  ///  @param Pos array containing co-ordinates of all the particles
  ///  @param mesh mesh related parameters -- connections and neighbours
  /// information
  /// @param lij_t0 initial distance between points of membrane
  /// @param is_attractive true if the zz sees the  bottom wall
  ///  @param mbrane  Membrane related parameters;
  /// @param mcpara Monte-Carlo related parameters
  /// @param AFM afm related parameter
  /// @return number of accepted moves

  int i, j, move;
  int nnbr_del1;
  int cm_idx_del1, cm_idx_del2;
  int cm_idx_add1, cm_idx_add2;
  int idx_del1, idx_del2;
  int idx_add1, idx_add2;
  int nframe;

  int nbr_add_1[12], nbr_add_2[12];
  int nbr_del_1[12], nbr_del_2[12];

  int N_nbr_del2, N_nbr_del1;
  int N_nbr_add2, N_nbr_add1;
  double det1, det2;
  Vec3d bef_ij, aft_ij;

  double KAPPA;
  bool yes, logic;

  nframe = get_nstart(mbrane.N, mbrane.bdry_type);
  move = 0;

  int idxn, up, down;

  for (i = 0; i < mcpara.one_mc_iter; i++) {
    // identify the pair to be divorced
    // stored as idx_del1 and idx_del2
    logic = false;
    while (!logic) {
      idx_del1 = RandomGenerator::intUniform(nframe, mbrane.N-1 );
      cm_idx_del1 = mesh.nghst * idx_del1;
      nnbr_del1 = mesh.numnbr[idx_del1];
      idxn = RandomGenerator::intUniform(0, mesh.nghst-1 );
      if (mesh.node_nbr_list[cm_idx_del1 + idxn] != -1) {
        idx_del2 = mesh.node_nbr_list[cm_idx_del1 + idxn];
        cm_idx_del2 = mesh.nghst * idx_del2;
        up = (idxn + 1 + nnbr_del1) % nnbr_del1;
        down = (idxn - 1 + nnbr_del1) % nnbr_del1;
        idx_add1 = mesh.node_nbr_list[idx_del1 * mesh.nghst + up];
        idx_add2 = mesh.node_nbr_list[idx_del1 * mesh.nghst + down];
        logic = idx_del2 > nframe && idx_add1 > nframe && idx_add2 > nframe;
      } else {
        logic = false;
      }
    }

    /* det1 = determinant(pos[idx_add1], pos[idx_add2], pos[idx_del1], mbrane.len); */
    /* det2 = determinant(pos[idx_add1], pos[idx_add2], pos[idx_del2], mbrane.len); */
    cm_idx_add1 = mesh.nghst * idx_add1;
    cm_idx_add2 = mesh.nghst * idx_add2;

        /* if (det1 * det2 < 0.0) { */
      aft_ij = pos[idx_add2] - pos[idx_add1];
      double dl = norm(aft_ij);
      N_nbr_del1 = mesh.numnbr[idx_del1];
      N_nbr_del2 = mesh.numnbr[idx_del2];
      N_nbr_add1 = mesh.numnbr[idx_add1];
      N_nbr_add2 = mesh.numnbr[idx_add2];
      bool flip_condt1, flip_condt2, flip_condt3;
      bool accept_flip;

      flip_condt1 = (dl < fl_para.fac_len_vertices*mbrane.av_bond_len);
      flip_condt2 =  N_nbr_del1 > fl_para.min_allowed_nbr && N_nbr_del2 > fl_para.min_allowed_nbr;
      flip_condt3 =  N_nbr_add1 < 9 && N_nbr_add2 < 9;

      accept_flip = flip_condt1 && flip_condt2 && flip_condt3;

      if (accept_flip) {
        move = move + 1;
        /* print_sanity(pos, mesh.node_nbr_list+cm_idx_del1,
         * mesh.node_nbr_list+cm_idx_del2, */
        /*         mesh.node_nbr_list+cm_idx_add1,
         * mesh.node_nbr_list+cm_idx_add2, */
        /*         idx_del1, idx_del2, idx_add1, idx_add2, (char*)"bef", i); */

        memcpy(nbr_del_1, &mesh.node_nbr_list[cm_idx_del1],
               sizeof(int) * mesh.nghst);
        memcpy(nbr_del_2, &mesh.node_nbr_list[cm_idx_del2],
               sizeof(int) * mesh.nghst);
        memcpy(nbr_add_1, &mesh.node_nbr_list[cm_idx_add1],
               sizeof(int) * mesh.nghst);
        memcpy(nbr_add_2, &mesh.node_nbr_list[cm_idx_add2],
               sizeof(int) * mesh.nghst);

        // form the bond
        N_nbr_add1 = add_nbr(nbr_add_1, mesh.numnbr[idx_add1], idx_add2,
                             idx_del1, idx_del2);
        N_nbr_add2 = add_nbr(nbr_add_2, mesh.numnbr[idx_add2], idx_add1,
                             idx_del1, idx_del2);

        // get divorced
        N_nbr_del1 = del_nbr(nbr_del_1, mesh.numnbr[idx_del1], idx_del2);
        N_nbr_del2 = del_nbr(nbr_del_2, mesh.numnbr[idx_del2], idx_del1);

        memcpy(mesh.node_nbr_list + cm_idx_del1, &nbr_del_1,
               sizeof(int) * mesh.nghst);
        memcpy(mesh.node_nbr_list + cm_idx_del2, &nbr_del_2,
               sizeof(int) * mesh.nghst);
        mesh.numnbr[idx_del1] = N_nbr_del1;
        mesh.numnbr[idx_del2] = N_nbr_del2;

        memcpy(mesh.node_nbr_list + cm_idx_add1, &nbr_add_1,
               sizeof(int) * mesh.nghst);
        memcpy(mesh.node_nbr_list + cm_idx_add2, &nbr_add_2,
               sizeof(int) * mesh.nghst);
        mesh.numnbr[idx_add1] = N_nbr_add1;
        mesh.numnbr[idx_add2] = N_nbr_add2;

            /* exit(0); */
      }
    /* } */
  }
  return move;
}
