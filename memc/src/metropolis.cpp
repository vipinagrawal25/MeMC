#include "metropolis.hpp"
#include "random_gen.hpp"
#include "hdf5_io.hpp"
#include <cmath>
#include <cstring>
#include <iomanip>
#include <mpi.h>

const double pi = 3.14159265358979323846264;
// #include <cstdio>
// #include <iomanip>
// #include <sstream>
// #include <unistd.h>

extern "C" void  MC_listread(char *, double *, double *, bool *,
                             int *, int *, bool *, bool *, int *, int *, double *, int *, char *);

// Remove this once fixed;
template<typename T>
string ZeroPadNumber(T num){
    ostringstream ss;
    ss << setw( 5 ) << setfill( '0' ) << (int)num;
    return ss.str();
}

int get_nstart(int, int);


int McP::initMC(int N, std::string fname){
  char tmp_fname[128], temp_algo[128];
  string parafile, outfile;

  parafile = fname+"/para_file.in";
  sprintf(tmp_fname, "%s", parafile.c_str());
  MC_listread(temp_algo, &dfac, &kBT, &is_restart,
              &tot_mc_iter, &dump_skip, &is_fluid, &is_semisolid, &min_allowed_nbr,
              &fluidize_every, &fac_len_vertices, &num_solid_points, tmp_fname);
  algo=temp_algo;
  one_mc_iter = 2*N;
  dfac = sqrt(8*pi/(2*N-4))/dfac;
  acceptedmoves = 0;

  solid_idx = (int *)calloc(N, sizeof(int));
  if(is_semisolid){
      if(num_solid_points > N){
          fprintf(stderr, "ERROR: num_solid_points (%d) exceeds total particles N (%d)\n",
                  num_solid_points, N);
          MPI_Abort(MPI_COMM_WORLD, 1);
      }
      int *index_solid = (int *)calloc(N, sizeof(int));
      string solidfile = fname + "/solid_index.h5";
      hdf5_io_read_int(index_solid, solidfile, "solid_idx");
      for(int i = 0; i < num_solid_points; i++)
          solid_idx[index_solid[i]] = 1;
      free(index_solid);
  }

  ofstream out_;
  out_.open( fname+"/mcpara.out");
  out_<< "# =========== monte carlo parameters ==========" << endl
      << " N " << N << endl
      << " algo = " << algo << endl
      << " dfac " << dfac << endl
      << " kbT " << kBT << endl
      << " is_restart " << is_restart << endl
      << " is_fluid " << is_fluid << endl
      << " is_semisolid " << is_semisolid << endl
      << " num_solid_points " << num_solid_points << endl
      << " tot_mc_iter " << tot_mc_iter << endl
      << " dump_skip " << dump_skip << endl
      << " min_allowed_nbr " << min_allowed_nbr << endl
      << " fluidize_every " << fluidize_every << endl;
  out_.close();

  return 1;
}
void McP::setEneVol(double radius) {
  EneMonitored = totEner;
  VolMonitored = totvol;
  volt0 = (4./3.)*pi*pow(radius,3);
}

double McP::evalEnergy(Vec3d *Pos, MESH_p mesh, std::fstream &fileptr, int itr){
  // BE bendeobj;
  // STE stretcheobj;
  double bende, stretche, pre=0;
  double stickener;
  double areat;
  double eself;
  //BDE 
  double bdryener=0.0;
if (fileptr.is_open()) {
  fileptr << itr << "  " << (double)acceptedmoves/(double)one_mc_iter<< "  "; 
  bende = beobj.bending_energy_total(Pos, mesh);
  stretche = steobj.stretch_energy_total(Pos, mesh);
  stickener = stickobj.stick_energy_total(Pos, mesh.N);

  if (bdeobj.do_bdry){
    bdryener=bdeobj.bde_total(Pos,mesh);
  }
  fileptr << bende << "  "<<stretche << "  "<< stickener << "  " << bdryener << " ";
  totvol = steobj.volume_total(Pos, mesh);
  // totarea = steobj.area_total(Pos, mesh);
  }

 totEner = bende+stretche+stickener+bdryener;

 if (steobj.dopressure()) {
   pre = steobj.PressureEnergyTotal(volt0, totvol);
   fileptr << pre << "  ";
   totEner += pre;
 }
 if (celllistobj.isSelfRepulsive()){
   eself = celllistobj.totalRepulsiveEnergy(Pos, mesh);
   fileptr << eself << "  ";
   totEner += eself;
 }
 areat = steobj.area_total(Pos, mesh);
 fileptr << totEner  << "  " << totvol  << "  " << areat  << "  " << volt0 << endl;

 if (is_fluid) {
   EneMonitored = totEner;
   VolMonitored = totvol;
 }
  return totEner;
}

// double McP::getarea(){return totarea;}
double McP::getvolume(){return totvol;}
bool McP::isrestart(){return is_restart;}
bool McP::isfluid(){return is_fluid;}
bool McP::issemisolid(){return is_semisolid;}
int* McP::getsolidIdx(){return solid_idx;}
void McP::mark_solid_neighbours(MESH_p mesh){
    for(int i = 0; i < mesh.N; i++){
        if(!solid_idx[i]) continue;
        int cm_idx = mesh.nghst * i;
        int num_nbr = mesh.numnbr[i];
        for(int j = 0; j < num_nbr; j++){
            int nbr = mesh.node_nbr_list[cm_idx + j];
            if(nbr >= 0) solid_idx[nbr] = 1;
        }
    }
}
int McP::fluidizeevery(){return fluidize_every;}
int McP::dumpskip(){return dump_skip;}
int McP::totaliter(){return tot_mc_iter;}
int McP::onemciter(){return one_mc_iter;}

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

bool McP::Boltzman(double DE, double activity) {
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
    yes = rand < exp(-DE / kBT);
  }
  return yes;
}
bool McP::Glauber(double DE, double activity) {
  /// @brief Glauber algorithm
  /// @param DE change in energy
  /// @param kbt boltzmann constant times temperature
  bool yes;
  double rand;
  DE += activity;
  rand = RandomGenerator::generateUniform(0.0,1.0);
  yes = rand < 1 / (1 + exp(DE / kBT));
  return yes;
}


double McP::energy_mc_3d(Vec3d *pos, MESH_p mesh, int idx) {
  double E_b, E_s, E_stick, E_afm, E_spr;
  int cm_idx, num_nbr;
  double Eself = 0.0;

  E_b = 0.0;
  E_s = 0.0;
  E_stick = 0.0;
  E_afm = 0.0;
  E_spr = 0.0;

  cm_idx = mesh.nghst * idx;
  num_nbr = mesh.numnbr[idx];

  E_b = beobj.bending_energy_ipart(pos, (int *)(mesh.node_nbr_list + cm_idx), num_nbr, idx);

  E_b += beobj.bending_energy_ipart_neighbour(pos, mesh, idx);

  E_s = steobj.stretch_energy_ipart(pos, (int *)(mesh.node_nbr_list + cm_idx),
                              num_nbr, idx, mesh.nghst);
  E_stick = stickobj.stick_energy_ipart(pos[idx], idx); 
  if(celllistobj.isSelfRepulsive()) Eself = celllistobj.computeSelfRep(pos, mesh, idx); 
  // cout << Eself <<endl;
//   if(st_p.do_stick)
//     if(afm.do_afm) E_afm = lj_afm(pos[idx], afm);

//     if(spring.do_spring) E_spr = spring_energy(pos[idx], idx, mesh, spring);
  return E_b + E_s + E_stick + E_afm + E_spr + Eself;
}
// //


int McP::monte_carlo_3d(Vec3d *pos, MESH_p mesh) {
  int i, num_nbr, cm_idx;
  double x_o, y_o, z_o, x_n, y_n, z_n;
  double de, Eini, Efin;
  double dxinc, dyinc, dzinc;
  double vol_i, vol_f;
  double dvol, de_vol, ini_vol, de_pressure;
  bool yes;
  int nframe;
  //
  nframe = get_nstart(mesh.N, mesh.bdry_type);
  acceptedmoves = 0;

  double bend_new[13]; // max nghst + 1

  for (i = 0; i < one_mc_iter; i++) {
    int idx = RandomGenerator::intUniform(nframe, mesh.N-1);
    cm_idx = idx*mesh.nghst;
    num_nbr = mesh.numnbr[idx];
    int *nbr_list = mesh.node_nbr_list + cm_idx;

    // --- Eini: use cached per-node bending values (no evaluation) ---
    double Eini_bend = beobj.bend_cache[idx];
    for(int k = 0; k < num_nbr; k++) Eini_bend += beobj.bend_cache[nbr_list[k]];
    double Eini_rest = steobj.stretch_energy_ipart(pos, nbr_list, num_nbr, idx, mesh.nghst)
                     + stickobj.stick_energy_ipart(pos[idx], idx);
    if(celllistobj.isSelfRepulsive()) Eini_rest += celllistobj.computeSelfRep(pos, mesh, idx);
    Eini = Eini_bend + Eini_rest;

    vol_i = steobj.volume_ipart(pos, nbr_list, num_nbr, idx);
    //
    x_o = pos[idx].x; y_o = pos[idx].y; z_o = pos[idx].z;
    //
    dxinc = (dfac) * (RandomGenerator::generateUniform(-1.0,1.0));
    dyinc = (dfac) * (RandomGenerator::generateUniform(-1.0,1.0));
    dzinc = (dfac) * (RandomGenerator::generateUniform(-1.0,1.0));
    //
    x_n = x_o + dxinc; y_n = y_o + dyinc; z_n = z_o + dzinc;
    //
    pos[idx].x = x_n; pos[idx].y = y_n; pos[idx].z = z_n;

    // --- Efin: compute fresh bending for idx and all its neighbours ---
    bend_new[0] = beobj.bending_energy_ipart(pos, nbr_list, num_nbr, idx);
    double Efin_bend = bend_new[0];
    for(int k = 0; k < num_nbr; k++){
        int nbr = nbr_list[k];
        int nbr_cm = nbr * mesh.nghst;
        bend_new[k+1] = beobj.bending_energy_ipart(pos,
                (int *)(mesh.node_nbr_list + nbr_cm), mesh.numnbr[nbr], nbr);
        Efin_bend += bend_new[k+1];
    }
    double Efin_rest = steobj.stretch_energy_ipart(pos, nbr_list, num_nbr, idx, mesh.nghst)
                     + stickobj.stick_energy_ipart(pos[idx], idx);
    if(celllistobj.isSelfRepulsive()) Efin_rest += celllistobj.computeSelfRep(pos, mesh, idx);
    Efin = Efin_bend + Efin_rest;

    de = (Efin - Eini);

    vol_f = steobj.volume_ipart(pos,
            (int *) (mesh.node_nbr_list + cm_idx), num_nbr, idx);
    dvol=0.5*(vol_f - vol_i);
    // std::cout << i << std::endl;
    if(steobj.dovol()){
    //   de_vol = vol_energy_change(mbrane, vol_p, dvol);
    //   // cout << de << "\t";
    //   de = (Efin - Eini) + de_vol;
    //     // cout << de << endl;
    }
    if(steobj.dopressure()){
      de_pressure = steobj.PV_change(dvol, volt0, VolMonitored);
      //de_pressure = steobj.getpressure()*log((VolMonitored + 2*dvol)/Volref);
      de = (Efin - Eini) + de_pressure;
    }
    double act = 0.0;
    if(actobj.is_active()){
      act = actobj.getActivityIdx(idx);
    }
    if (algo == "mpolis"){
      yes = Boltzman(de, act);
    } else if (algo == "glauber") {
      yes = Glauber(de, act);
    }
    //
    if(yes) {
      acceptedmoves +=  1;
      EneMonitored += de;
      VolMonitored += 2*dvol;
      // update bend cache for idx and all its neighbours
      beobj.bend_cache[idx] = bend_new[0];
      for(int k = 0; k < num_nbr; k++)
        beobj.bend_cache[nbr_list[k]] = bend_new[k+1];
    } else {
      pos[idx].x = x_o;
      pos[idx].y = y_o;
      pos[idx].z = z_o;
    }
  }
  return acceptedmoves;
}
// //

int McP::monte_carlo_fluid(Vec3d *pos, MESH_p mesh, double av_bond_len) {

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
  bool yes, logic;

  nframe = get_nstart(mesh.N, mesh.bdry_type);
  move = 0;

  int idxn, up, down;

  for (i = 0; i < one_mc_iter; i++) {
    // identify the pair to be divorced
    // stored as idx_del1 and idx_del2
    logic = false;
    while (!logic) {
      idx_del1 = RandomGenerator::intUniform(nframe, mesh.N-1 );
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
        if(is_semisolid && logic)
            logic = solid_idx[idx_del1] + solid_idx[idx_del2] +
                    solid_idx[idx_add1] + solid_idx[idx_add2] == 0;
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

      flip_condt1 = (dl < fac_len_vertices*av_bond_len);
      flip_condt2 =  N_nbr_del1 > min_allowed_nbr && N_nbr_del2 > min_allowed_nbr;
      flip_condt3 =  N_nbr_add1 < 9 && N_nbr_add2 < 9;

      accept_flip = flip_condt1 && flip_condt2 && flip_condt3;
      if (accept_flip) {
        move = move + 1;

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
      }
  }
  return move;
}

//adding boundary sampling loop 
//for stress-free boundary
int McP::monte_carlo_bdry(Vec3d *pos, MESH_p mesh){
  if (!bdeobj.do_bdry) return 0; //if boundary is not stress-free
  int num_nbr,cm_idx;
  double x_o,y_o,z_o,x_n,y_n,z_n;
  double de,Eini,Efin;
  double dxinc,dyinc,dzinc;
  bool yes;
  int bdry_acceptedmoves=0;

  int nframe=bdeobj.num_bdry_nodes;
  double bend_new[13];

  double dfac_bdry=dfac; //is this supposed to be the same?

  for (int i=0;i<2*nframe;i++){  //since bulk mc is 2*N per mciter 

    int idx=RandomGenerator::intUniform(0,nframe-1);

    if (bdeobj.is_clamped_vertex(idx)) continue;

    cm_idx=idx*mesh.nghst;
    num_nbr=mesh.numnbr[idx];
    int *nbr_list=mesh.node_nbr_list+cm_idx;

    //bending energy change for neighbours (non-boundary)
    double Eini_bend=beobj.bend_cache[idx];
    for (int k=0; k<num_nbr; k++){
       Eini_bend+=beobj.bend_cache[nbr_list[k]];
    }
    double Eini_rest=steobj.stretch_energy_ipart(pos,nbr_list, num_nbr, idx,mesh.nghst)+stickobj.stick_energy_ipart(pos[idx],idx);
    if (celllistobj.isSelfRepulsive()){
      Eini_rest+=celllistobj.computeSelfRep(pos,mesh,idx);
    }
    
    double Eini_gc=bdeobj.bde_ipart(pos,mesh,idx); //change due to geodesic curvature integral
    Eini=Eini_bend+Eini_rest+Eini_gc;

    //propose update
    x_o=pos[idx].x;
    y_o=pos[idx].y;
    z_o=pos[idx].z;
  
    dxinc=dfac_bdry*RandomGenerator::generateUniform(-1.0,1.0);
    dyinc=dfac_bdry*RandomGenerator::generateUniform(-1.0,1.0);
    dzinc=dfac_bdry*RandomGenerator::generateUniform(-1.0,1.0);
    
    x_n = x_o + dxinc; y_n = y_o + dyinc; z_n = z_o + dzinc;
    //
    pos[idx].x = x_n; pos[idx].y = y_n; pos[idx].z = z_n;

    bend_new[0]=beobj.bending_energy_ipart(pos,nbr_list,num_nbr,idx);
    double Efin_bend=bend_new[0];

    for (int k=0; k<num_nbr; k++){
      int nbr=nbr_list[k];
      int nbr_cm=nbr*mesh.nghst;
      bend_new[k+1]=beobj.bending_energy_ipart(pos,(int*)(mesh.node_nbr_list+nbr_cm),mesh.numnbr[nbr],nbr);
      Efin_bend+=bend_new[k+1];
    }

    double Efin_rest=steobj.stretch_energy_ipart(pos,nbr_list,num_nbr,idx,mesh.nghst)+stickobj.stick_energy_ipart(pos[idx],idx);
    if (celllistobj.isSelfRepulsive()){
      Efin_rest+=celllistobj.computeSelfRep(pos,mesh,idx);
    }
   
    double Efin_gc=bdeobj.bde_ipart(pos,mesh,idx);
    Efin=Efin_bend+Efin_rest+Efin_gc;

    de=(Efin-Eini);


    double act=0.0;
    if (actobj.is_active()){
       act=actobj.getActivityIdx(idx);
    }

    if (algo=="mpolis"){
      yes=Boltzman(de,act);
    }
    else if (algo=="glauber"){
      yes=Glauber(de,act);
    }


    if (yes){
      bdry_acceptedmoves+=1;
      EneMonitored+=de;
      beobj.bend_cache[idx]=bend_new[0];
      for (int k=0; k<num_nbr; k++){
        beobj.bend_cache[nbr_list[k]]=bend_new[k+1];
      }
    } else {
      pos[idx].x=x_o;
      pos[idx].y=y_o;
      pos[idx].z=z_o;
    }
  }
  return bdry_acceptedmoves;
}


    
