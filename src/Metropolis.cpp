#include "Vector.h"
#include "global.h"
#include "subroutine.h"
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <sstream>
#include <unistd.h>

/* void init_rng(uint32_t seed_val) { */

/*   ///  @brief Generates random number */
/*   rng.seed(seed_val); */
/* } */

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
  /* std::uniform_real_distribution<> rand_real(0, 1); */
  yes = (DE <= 0.e0);
  if (!yes) {
    rand = RandomGenerator::generateUniform(0.0, 1.0); //  rand_real(rng);
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
  rand = RandomGenerator::generateUniform(0.0, 1.0); //  rand_real(rng);
  yes = rand < 1 / (1 + exp(DE / mcpara.kBT));
  return yes;
}

/* RandomGenerator::generateUniform */
double rand_inc_theta(double th0, double dfac) {
  /// @brief increment the polar angle randomly
  double dth;
  double tmp_th0;
 // std::uniform_real_distribution<> rand_real(-1, 1);

  tmp_th0 = 10;
  while (tmp_th0 > pi || tmp_th0 < 0) {
    dth = (pi / dfac) * RandomGenerator::generateUniform(0.0,1.0); 
    //(rand_real(rng));
    tmp_th0 = th0 + dth;
  }
  return dth;
}

double energy_mc_3d(Vec3d *pos, Vec3d *pos_t0, MESH_p mesh, double *lij_t0, 
                    int idx, double *area_i, MBRANE_p mbrane, AREA_p area_p, STICK_p st_p, 
                    VOL_p vol_p, AFM_p afm) {
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
  double area_idx;
  Vec2d be_ar;
  int cm_idx, num_nbr;

  int nframe;
  nframe = get_nstart(mbrane.N, mbrane.bdry_type);

  E_b = 0.0;
  E_s = 0.0;
  E_stick = 0.0;
  E_afm = 0.0;

  cm_idx = mesh.nghst * idx;
  num_nbr = mesh.numnbr[idx];

  be_ar = bending_energy_ipart(pos, (int *)(mesh.node_nbr_list + cm_idx), num_nbr,
                             idx, mbrane);
  E_b = be_ar.x;
  /* area_idx = be_ar.y; */

  be_ar = bending_energy_ipart_neighbour(pos, mesh, idx, mbrane);

  E_b += be_ar.x;

  if(area_p.is_tether){
      E_s = stretch_energy_ipart(pos, (int *)(mesh.node_nbr_list + cm_idx),
              (lij_t0 + cm_idx),  num_nbr, idx, area_p);
  }

  /* *area_i = area_idx + be_ar.y; */ 

  area_idx = area_ipart(pos,  (int *) (mesh.node_nbr_list + cm_idx),
                 num_nbr, idx);
   *area_i = area_idx;


      if(st_p.do_stick && st_p.is_attractive[idx])
          E_stick = stick_bottom_surface(pos[idx], pos_t0[idx], st_p); 

      if(afm.do_afm) E_afm = lj_afm(pos[idx], afm);
      /* fprintf(stderr, "%d \n", idx); */

  return E_b + E_s + E_stick + E_afm;
}

int monte_carlo_3d(Vec3d *pos, Vec3d *pos_t0, MESH_p mesh, double *lij_t0, 
                   MBRANE_p mbrane, MC_p mcpara, AREA_p area_p,  STICK_p st_p,
                   VOL_p vol_p, AFM_p afm,
                   ACTIVE_p activity, SHEAR_p shear) {

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
  double kappa;
  double area_i, area_f, ini_ar;
  double d_ar, de_area;
  bool yes;

  int nframe = 0;
  nframe = get_nstart(mbrane.N, mbrane.bdry_type);
  /* std::uniform_int_distribution<uint32_t> rand_int(nframe, mbrane.N - 1); */
  /* std::uniform_real_distribution<> rand_real(-1, 1); */
  ini_vol = (4. / 3.) * pi * pow(mbrane.radius, 3);
  /* if(mbrane.bdry_type == 0 || mbrane.bdry_type == 1 ) */
  ini_ar = 4 * pi * mbrane.radius*mbrane.radius; 
  kappa = vol_p.coef_vol_expansion;
  move = 0;

  for (i = 0; i < mcpara.one_mc_iter; i++) {
     int idx =  RandomGenerator::intUniform(nframe, mbrane.N - 1);
      /* int idx = rand_int(rng); */
      cm_idx = mesh.nghst * idx;
      num_nbr = mesh.numnbr[idx];
      Eini = energy_mc_3d(pos, pos_t0, mesh, lij_t0,  idx, &area_i, mbrane, area_p, st_p, vol_p,
              afm);
      if(vol_p.do_volume) vol_i = volume_ipart(pos,
              (int *) (mesh.node_nbr_list + cm_idx), num_nbr, idx);

      x_o = pos[idx].x;
      y_o = pos[idx].y;
      z_o = pos[idx].z;

      dxinc = (mcpara.delta / mcpara.dfac) * (RandomGenerator::generateUniform(-1,1));
      dyinc = (mcpara.delta / mcpara.dfac) * (RandomGenerator::generateUniform(-1,1));
      dzinc = (mcpara.delta / mcpara.dfac) * (RandomGenerator::generateUniform(-1,1));

      if(!area_p.is_tether){
          double rand_random_ = (RandomGenerator::generateUniform(-1,1));
          double magr_ = norm(pos[idx]);
          dxinc = (mcpara.delta / mcpara.dfac) * rand_random_*pos[idx].x/magr_;
          dyinc = (mcpara.delta / mcpara.dfac) * rand_random_*pos[idx].y/magr_;
          dzinc = (mcpara.delta / mcpara.dfac) * rand_random_*pos[idx].z/magr_;
      }

      x_n = x_o + dxinc;
      y_n = y_o + dyinc;
      z_n = z_o + dzinc;

      pos[idx].x = x_n;
      pos[idx].y = y_n;
      pos[idx].z = z_n;

      Efin = energy_mc_3d(pos, pos_t0, mesh, lij_t0,  idx, &area_f, mbrane, area_p, st_p, vol_p,
              afm);

      de = (Efin - Eini);
      if(!area_p.is_tether){
          d_ar = area_f - area_i;
         de_area = (2*d_ar/(ini_ar*ini_ar))*(mbrane.area[0]  - ini_ar)
                  + (d_ar/ini_ar)*(d_ar/ini_ar);

            de +=  area_p.Ka*de_area;
      }
      if(vol_p.do_volume){
          vol_f = volume_ipart(pos,
                  (int *) (mesh.node_nbr_list + cm_idx), num_nbr, idx);
          dvol=(vol_f - vol_i);

          if(!vol_p.is_pressurized){
              de_vol = (2*dvol/(ini_vol*ini_vol))*(*mbrane.volume  - ini_vol)
                  + (dvol/ini_vol)*(dvol/ini_vol);
              de +=  kappa*de_vol;
          }
          if(vol_p.is_pressurized){
              de_pressure = vol_p.pressure*dvol;
              /* de_pressure = vol_p.pressure*dvol*(*mbrane.volume/ini_vol); */
              de +=  de_pressure;
          }
      }
     /* yes = Metropolis(de, 0.0, mcpara); */
      if (mcpara.algo == "mpolis") {
          yes = Metropolis(de, activity.activity[idx], mcpara);
      } else if (mcpara.algo == "glauber") {
          yes = Glauber(de, activity.activity[idx], mcpara);
      }

      if (yes) {
          move = move + 1;
          *mbrane.tot_energy += de;
          if(vol_p.do_volume) *mbrane.volume += dvol; 
          *mbrane.area += d_ar; 
      } else {
          pos[idx].x = x_o;
          pos[idx].y = y_o;
          pos[idx].z = z_o;
      }
  }
  return move;
}

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
    int idx = RandomGenerator::intUniform(n_ghost, para.N-1);
    Eini = pairlj_ipart_energy(Pos, neib[idx].list, neib[idx].cnt, idx, para,
                               metric);
    x_o = Pos[idx].x;
    y_o = Pos[idx].y;

    if (is_cart) {
      dxinc = (para.sigma / mcpara.dfac) * RandomGenerator::generateUniform(-1,1); 
      x_n = fmod((x_o + dxinc + 30 * para.len), para.len);
      Pos[idx].x = x_n;

      dyinc = (para.sigma / mcpara.dfac) * RandomGenerator::generateUniform(-1,1);
      y_n = fmod((y_o + dyinc + 30 * para.len), para.len);
      Pos[idx].y = y_n;
    }
    if (is_sph) {
      dxinc = rand_inc_theta(Pos[idx].x, mcpara.dfac);
      dyinc = (para.sigma / mcpara.dfac) * RandomGenerator::generateUniform(-1,1);
      x_n = x_o + para.sigma * dxinc;
      y_n = fmod((y_o + dyinc + 30 * 2 * pi), 2 * pi);
      Pos[idx].x = x_n;
      Pos[idx].y = y_n;
    }

    Efin = pairlj_ipart_energy(Pos, neib[idx].list, neib[idx].cnt, idx, para,
                               metric);
    de = (Efin - Eini);
    if (Metropolis(de, 0.0, mcpara)) {
      move = move + 1;
    } else {
      Pos[idx].x = x_o;
      Pos[idx].y = y_o;
    }
  }
  return move;
}

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

  double kappa;
  bool yes, logic;
  bool do_move_pt;

  nframe = get_nstart(mbrane.N, mbrane.bdry_type);

  /*   std::uniform_int_distribution<uint32_t> rand_int(nframe, mbrane.N - 1); */
  /*   std::uniform_int_distribution<uint32_t> rand_nbr(0, mesh.nghst - 1); */
  /*   std::uniform_real_distribution<> rand_real(-1, 1); */

  move = 0;

  int idxn, up, down;

  for (i = 0; i < mcpara.one_mc_iter; i++) {
    // identify the pair to be divorced
    // stored as idx_del1 and idx_del2
    logic = false;
    while (!logic) {
      idx_del1 = RandomGenerator::intUniform(nframe, mbrane.N-1); //rand_int(rng);
      cm_idx_del1 = mesh.nghst * idx_del1;
      nnbr_del1 = mesh.numnbr[idx_del1];
      idxn = RandomGenerator::intUniform(0, mesh.nghst-1); //rand_nbr(rng);
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

    if(fl_para.is_semifluid){
      do_move_pt = fl_para.solid_idx[idx_del1] + fl_para.solid_idx[idx_del2] + 
        fl_para.solid_idx[idx_add1] + fl_para.solid_idx[idx_add2] == 0;
    }else{
      do_move_pt = true;
    }
    /* det1 = determinant(pos[idx_add1], pos[idx_add2], pos[idx_del1], mbrane.len); */
    /* det2 = determinant(pos[idx_add1], pos[idx_add2], pos[idx_del2], mbrane.len); */
    if(do_move_pt){

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

      /* cout<< N_nbr_add2 << "  " << flip_condt2<< "  " << N_nbr_add1 << endl; */
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
    }
  }
  return move;
}
