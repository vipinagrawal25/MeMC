#include "metropolis.hpp"
#include "random_gen.hpp"
#include "multicomp.hpp"
#include "electrostatics.hpp"
#include "selfavoidance.hpp"
#include "linetension.hpp"

#include <cmath>
#include <cstring>
#include <tuple>
#include <functional>
#include <algorithm>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <mpi.h>

const double pi = 3.14159265358979323846264;

using namespace std;
extern "C" void  MC_listread(char *, double *, double *, bool *, int*, int *, 
   int *, bool *, int *, int *, double *, bool *, int*, char *);
int get_nstart(int, int);

McP::McP (BE &beobj, STE &steobj, MulCom &lipidobj, ESP &chargeobj,
   SelfAvoid &repulsiveobj, LTN &lineobj):
beobj(beobj), steobj(steobj), lipidobj(lipidobj), chargeobj(chargeobj), 
repulsiveobj(repulsiveobj), lineobj(lineobj){};
//
// void McP::startcycle(int cycle){
//    if (0 <= cycle && cycle < 3){
//       this->kBT = initial_kBT;
//       this->dfac = initial_dfac;
//    }else{
//       this->kBT = initial_kBT * pow(10.0, -(cycle - 4));
//       this->dfac = initial_dfac / pow(2.0, (cycle - 4));
//    }
//    cout << "Cycle: " << cycle << " kBT: " << kBT << " dfac: " << dfac << endl;
// }

void McP::startcycle(int cycle){

    // First two cycles start at initial temperature
    if (cycle < 2) {
        kBT = initial_kBT;   // 0.1
    }
    else {
        // cycle 2 → 1e-2
        // cycle 3 → 1e-3
        kBT = initial_kBT * pow(10.0, -(cycle - 1));
    }

    dfac = initial_dfac;

    cout << "Cycle: " << cycle 
         << " starting kBT: " << kBT 
         << endl;
}

//
void McP::updateparam(int anneal, string fname){
   if (anneal>0){
      dfac=dfac/2;
      kBT=kBT*0.1;

   ofstream out_;
   out_.open( fname+"/mcpara.out", ios::app );
   out_<< "# =========== parameter after annealing = " << anneal <<  " ==========" << endl
      << " dfac " << dfac << endl
      << " kbT " << kBT << endl
      << " is_evalEnergyrestart " << is_restart << endl
      << " tot_mc_iter " << tot_mc_iter << endl
      << " dump_skip " << dump_skip << endl;
   out_.close();
   }
}
//
int McP::initMC(MESH_p mesh, string fname){
   int N = mesh.N;
   double radius = mesh.radius;
   // string topology=mesh.topology;
   char tmp_fname[128];
   char temp_algo[128] = {0};
   string parafile, outfile;
   parafile = fname+"/para_file.in";
   sprintf(tmp_fname, "%s", parafile.c_str());
   MC_listread(temp_algo, &dfac, &kBT, &is_restart, &nanneal_cycle,
              &tot_mc_iter, &dump_skip, &is_fluid, &min_allowed_nbr,
              &fluidize_every, &fac_len_vertices, &iexch, &nexch_iter, tmp_fname);

   std::string algo_str(temp_algo, 128);
   size_t nullpos = algo_str.find('\0');
   if (nullpos != std::string::npos) {
      algo_str.resize(nullpos);
   }
   auto trim_whitespace = [](std::string &s){
      const auto ws = " \t\r\n";
      size_t first = s.find_first_not_of(ws);
      if (first == std::string::npos) {
         s.clear();
         return;
      }
      size_t last = s.find_last_not_of(ws);
      s = s.substr(first, last - first + 1);
   };
   trim_whitespace(algo_str);

   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   // std::cout << "RANK " << rank << " Algo = [" << algo_str << "] len=" << algo_str.size() << std::endl;

   if (algo_str.empty()) {
      throw std::runtime_error("Empty algo string on MPI rank " + std::to_string(rank));
   }

   // strict validation
   if (algo_str != "mpolis" && algo_str != "Glauber") {
      throw std::runtime_error("Invalid algo string on rank " + std::to_string(rank) + ": [" + algo_str + "]");
   }

   ini_tot_mc_iter = tot_mc_iter;
   one_mc_iter = 2*N;
   dfac=mesh.av_bond_len/dfac;
   initial_dfac = dfac;
   initial_kBT = kBT;
   repulsiveobj.setKBT(kBT);
   acceptedmoves = 0;
   if (mesh.ncomp==1) iexch=false;
   if (!chargeobj.isexch() && !beobj.isexch() && !lineobj.calculate()) iexch=false;
   ofstream out_;
   out_.open( fname+"/mcpara.out" );
   out_ << "# =========== monte carlo parameters ==========" << endl
        << " N " << N << endl
        << " dfac " << dfac << endl
        << " kbT " << kBT << endl
        << " is_restart " << is_restart << endl
        << " is_fluid " << is_fluid << endl
        << " nanneal cycle = " << nanneal_cycle << endl
        << " tot_mc_iter " << tot_mc_iter << endl
        << " dump_skip " << dump_skip << endl
        << " min_allowed_nbr " << min_allowed_nbr << endl
        << " fluidize_every " << fluidize_every << endl
        << " iexch " << iexch << endl
        << " Number_exch_iter " << nexch_iter * one_mc_iter << endl;

   if (chargeobj.calculate() && lineobj.calculate()){
      out_ << " Energy_mc_3d = energy_mc_bestchli "<< endl;
      energy_mc_3d = [this](vector<double>& vec, Vec3d* vec_ptr, MESH_p mesh,
                     int val, int val2, int val3, BoundaryType btype ) -> double {
      return this->energy_mc_bestchli(vec, vec_ptr, mesh, val, val2, val3, btype);};
   }else if (chargeobj.calculate() && repulsiveobj.isSelfRepulsive()){
      out_ << " Energy_mc_3d = energy_mc_bestchrep "<< endl;
      energy_mc_3d = [this](vector<double>& vec, Vec3d* vec_ptr, MESH_p mesh, 
                  int val, int val2, int val3, BoundaryType btype ) -> double {
      return this->energy_mc_bestchrep(vec, vec_ptr, mesh, val, val2, val3, btype);};
   }else if(chargeobj.calculate()){
      out_ << " Energy_mc_3d = energy_mc_bestch" << endl;
      energy_mc_3d = [this](vector<double>& vec, Vec3d* vec_ptr, MESH_p mesh, int val,
                           int val2, int val3, BoundaryType btype) -> double {
      return this->energy_mc_bestch(vec , vec_ptr, mesh, val, val2, val3, btype);};
   }else{
      out_ << " Energy_mc_3d = energy_mc_best" << endl;
      energy_mc_3d = [this](vector<double>& vec, Vec3d* vec_ptr, MESH_p mesh, 
                    int val, int val2, int val3, BoundaryType btype) -> double {
      return this->energy_mc_best(vec, vec_ptr, mesh, val, val2, val3, btype);};
   }
   
   algo = algo_str;
   if (algo == "mpolis"){
      out_ << " Algo = Metropolis"<< endl;
      Algo = [this](double DE, double activity) -> double {
      return this->Boltzman(DE, activity);};
   }else if(algo == "Glauber"){
      out_ << " Algo = Glauber"<< endl;
      Algo = [this](double DE, double activity) -> double {
      return this->Glauber(DE, activity);};
   }else {
      throw std::runtime_error(
         "Unknown algorithm AFTER parsing: [" + algo + "]"
      );
   }
   // std::cout << "INIT OK: Algo assigned = " << algo << std::endl;

   if (!Algo) {
      throw std::runtime_error("FATAL: Algo not assigned after parsing");
   }

   if (lineobj.calculate()){
      out_ << " Energy_mc_exch = energy_mc_bechli" << endl;
      energy_mc_exch = [this](vector<double>& vec, Vec3d* vec_ptr, MESH_p mesh,
                  int val, int val2,
                  int val3, int val4, int val5, int val6,
                  BoundaryType btype1, BoundaryType btype2) -> double {
      return this->energy_mc_bechli(vec, vec_ptr, mesh, val, val2, val3,
                                 val4, val5, val6, btype1, btype2);};
   }else if(chargeobj.isexch() && beobj.isexch()){
      out_ << " Energy_mc_exch = energy_mc_bech" << endl;
      energy_mc_exch = [this](vector<double>& vec, Vec3d* vec_ptr, MESH_p mesh,
                  int val, int val2, 
                  int val3, int val4, 
                  int val5, int val6,
                  BoundaryType btype1, BoundaryType btype2) -> double {
      return this->energy_mc_bech(vec , vec_ptr, mesh, val, val2, val3, 
                                 val4, val5, val6, btype1, btype2);};
   }else if(chargeobj.isexch()){
      out_ << " Energy_mc_exch = energy_mc_ch" << endl;
      energy_mc_exch = [this](vector<double>& vec, Vec3d* vec_ptr, MESH_p mesh,
                  int val, int val2,
                  int val3, int val4,
                  int val5, int val6,
                  BoundaryType btype1, BoundaryType btype2) -> double {
      return this->energy_mc_ch(vec , vec_ptr, mesh, val, val2, val3, 
                                 val4, val5, val6);};
   }else if(beobj.isexch()){
      out_ << " Energy_mc_exch = energy_mc_be" << endl;
      energy_mc_exch = [this](vector<double>& vec, Vec3d* vec_ptr, MESH_p mesh,
                  int val, int val2,
                  int val3, int val4,
                  int val5, int val6,
                  BoundaryType btype1, BoundaryType btype2) -> double {
      return this->energy_mc_be(vec , vec_ptr, mesh, val, val2, val3, 
                  val4, val5, val6, btype1, btype2);};
   }

   if (iexch && !energy_mc_exch) {
      throw std::runtime_error("FATAL: energy_mc_exch not initialized in initMC() for exchange-enabled run");
   }
   
   if (exchtype=="Global" || exchtype == "global"){
      out_ << "Component exchange type = Global" << endl;
      get_idx2 = [this](int num_nbr, int *node_nbr_list, int cm_idx, int nframe,
         int N) -> int {return this->global_idx(nframe, N);};
   }else{
      out_ << "Component exchange type = Local" << endl;
      get_idx2 = [this](int num_nbr, int *node_nbr_list, int cm_idx, int nframe,
         int N) -> int {return this->local_idx(num_nbr, node_nbr_list, cm_idx);};
   }

   out_.close();
   volt0=mesh.ini_vol;
   if (!Algo) {
      throw std::runtime_error("FATAL: Algo not initialized in initMC()");
   }
   if (iexch && !energy_mc_exch) {
      throw std::runtime_error("FATAL: energy_mc_exch not initialized in initMC()");
   }
   if (iexch && !get_idx2) {
      throw std::runtime_error("FATAL: get_idx2 not initialized in initMC()");
   }
   return 1;
}
/*-----------------------*/
void McP::setEneVol() {
   EneMonitored = totEner;
   VolMonitored = totvol;
}
/*-----------------------*/
double McP::evalEnergy(MESH_p mesh){
   // if (fileptr.is_open()) {
   // fileptr << itr << "  " << (double)acceptedmoves/(double)one_mc_iter<< "  ";
   Vec3d *Pos = mesh.pos;
   bende = beobj.bending_energy_total(Pos, mesh);
   stretche = steobj.stretch_energy_total(Pos, mesh);
   
   	if (steobj.doarea()){
      stretche =steobj.area_energy_total(mesh);
      totEner += stretche;
   	}
   	totEner = bende+stretche;
   	if (chargeobj.calculate()){
      electroe = chargeobj.debye_huckel_total(Pos, mesh.N);
      totEner += electroe;
   	}
   	// It's okay to have if statement here.
   	if (lipidobj.calculate()) {
      regsole = lipidobj.reg_soln_tot(Pos, mesh);
      totEner+=regsole;
   	}
   	if (mesh.sphere) totvol = steobj.volume_total(Pos, mesh);
   	if (steobj.dovol()) {
      vole = steobj.getkappa()*(VolMonitored/volt0-1)*(VolMonitored/volt0-1);
      totEner += vole;
   	}
   	if (steobj.dopressure()){
      pre = -steobj.getpressure() * totvol;
      totEner += pre;
   	}
   	if (repulsiveobj.isSelfRepulsive()){
      selfe = repulsiveobj.totalRepulsiveEnergy(mesh);
      totEner += selfe;
   }
   if (lineobj.calculate()>0){
      linee = lineobj.energy_total();
      totEner += linee;
   }
   EneMonitored = totEner;
   VolMonitored = totvol;
   // totarea = steobj.area_total(mesh);
   AreaMonitored = steobj.area_total(mesh);
   return totEner;
}
/*----------------------------------------------------------*/
void McP::wHeader(const MESH_p &mesh, std::fstream &fid){
    std::string log_headers = "#iter acceptedmoves kBT bend_e stretch_e ";
    if (chargeobj.calculate()) log_headers+="electroe ";
    if(steobj.dopressure()) {log_headers+=" Pressure_e ";}
    if(steobj.dovol()) {log_headers+=" Volume_e ";}
    if (repulsiveobj.isSelfRepulsive()) {log_headers+=" Repulsive_e ";}
    if (lineobj.calculate()>0){log_headers+=" Line_e ";}
    log_headers+="total_e ";
    if (mesh.sphere){log_headers+="Volume ";}
    log_headers+="area";
    fid << log_headers << endl;
}
/*-----------------------*/
void McP::write_energy(fstream &fileptr, int itr, const MESH_p &mesh){
   if (fileptr.is_open()){
      fileptr << itr << " " << (double)acceptedmoves/(double)one_mc_iter<< "  ";
      fileptr << kBT << " ";
      fileptr << bende << " " << stretche << "  ";
      if (chargeobj.calculate()) fileptr << electroe << " ";
      if (lipidobj.calculate()) fileptr << regsole << " ";
      if (steobj.dopressure()) fileptr << pre << "  ";
      if (steobj.dovol()) fileptr << vole << " ";
      if (repulsiveobj.isSelfRepulsive()) fileptr << selfe << " ";
      if (lineobj.calculate()) fileptr << linee << " ";
      fileptr << EneMonitored  << "  ";
      if(mesh.sphere) fileptr << VolMonitored << " ";
      fileptr << AreaMonitored  << endl;
   }
}
/*-----------------------*/
double McP::getarea(){return AreaMonitored;}
double McP::getvolume(){return totvol;}
bool McP::isrestart(){return is_restart;}
bool McP::isfluid(){return is_fluid;}
int McP::fluidizeevery(){return fluidize_every;}
int McP::dumpskip(){return dump_skip;}
int McP::totaliter(){return tot_mc_iter;}
int McP::onemciter(){return one_mc_iter;}
// int McP::monitoredVol(){return VolMonitored;}
// int McP::monitoredEn(){return EneMonitored;}
//
int del_nbr(int *nbrs, int numnbr, int idx){
  // delet int idx between i1 and i2 in the nbrs list
  int new_numnbr, delete_here;
  bool logic;

  logic = false;
  delete_here = 0;

  while (!logic) {
    logic = (nbrs[delete_here] == idx);
    ++delete_here;
  }

  memcpy(nbrs + delete_here - 1, &nbrs[delete_here],
         sizeof(int) * (numnbr - delete_here + 1));

  return numnbr - 1;
}

int add_nbr(int *nbrs, int numnbr, int idx, int i1, int i2){
  // add int idx between i1 and i2 in the nbrs list
   int insert_here;
   bool logic;

   logic = false;
   insert_here = 0;

   while (!logic) {
      logic = (nbrs[insert_here] == i1) || (nbrs[insert_here] == i2);
      ++insert_here;
   }

   logic = (nbrs[insert_here] == i1) || (nbrs[insert_here] == i2);
   if (logic) {
      memcpy(nbrs + insert_here, &nbrs[insert_here - 1],
           sizeof(int) * (numnbr - insert_here + 1));
      nbrs[insert_here] = idx;
   }else {
      insert_here = 0;
      memcpy(nbrs + insert_here, &nbrs[insert_here - 1],
           sizeof(int) * (numnbr - insert_here + 1));
      nbrs[insert_here] = idx;
   }

   return numnbr + 1;
}
//
bool McP::Boltzman(double DE, double activity){
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
   // if(yes)cout << DE << endl;
   if (!yes) {
      rand = RandomGenerator::generateUniform(0.0,1.0);
      yes = rand < exp(-DE / kBT);
   }
   return yes;
}
//
bool McP::Glauber(double DE, double activity){
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
//
inline double McP::energy_mc_best(vector<double> &energy, Vec3d *pos, MESH_p mesh,
				int idx, int cm_idx, int num_nbr, BoundaryType btype){
   int *nbrcm=mesh.node_nbr_list + cm_idx;
	double lijsq[num_nbr];
	energy[0]  = beobj.bending_energy_ipart(pos, nbrcm, num_nbr, idx, mesh.boxlen, btype,lijsq);
	energy[0] += beobj.bending_energy_ipart_neighbour(pos, mesh, idx);
	if (steobj.getyy1()!=0 && steobj.getyy2()!=0){
		energy[1] = steobj.stretch_energy_ipart(lijsq, num_nbr, idx, mesh.nghst);
	}
	if (steobj.doarea()){
		energy[1] =  steobj.area_energy_ipart(pos, nbrcm, num_nbr, idx, mesh.boxlen, btype);
	}
	return energy[0] + energy[1];
}
//
inline double McP::energy_mc_bestch(vector<double> &energy, Vec3d *pos, MESH_p mesh, int idx, int cm_idx, int num_nbr, BoundaryType btype){
   double Etot=energy_mc_best(energy, pos, mesh, idx, cm_idx, num_nbr, btype);
   energy[2] = chargeobj.debye_huckel_ipart(pos, idx, mesh.N);
   return Etot + energy[2];
}
//
inline double McP::energy_mc_bestchli(vector<double> &energy, Vec3d *pos, MESH_p mesh, 
         int idx, int cm_idx, int num_nbr, BoundaryType btype){
   double Etot=energy_mc_bestch(energy, pos, mesh, idx, cm_idx, num_nbr, btype);
   energy[4] = lineobj.energy_ipart(idx);
   return Etot + energy[4];
}
//
inline double McP::energy_mc_bestchrep(vector<double> &energy, Vec3d *pos, MESH_p mesh,
            int idx, int cm_idx, int num_nbr, BoundaryType btype){
   double Etot=energy_mc_bestch(energy, pos, mesh, idx, cm_idx, num_nbr, btype);
   energy[3] = repulsiveobj.computeSelfRep(mesh, idx);
   return Etot + energy[3];
}
//
inline double McP::energy_mc_bech(vector<double> &energy, Vec3d *pos, MESH_p mesh,
    	int idx1, int idx2, int cm_idx1, int cm_idx2, int num_nbr1, int num_nbr2,
		BoundaryType btype1, BoundaryType btype2){
   int *nbrcm1=mesh.node_nbr_list + cm_idx1;
   int *nbrcm2=mesh.node_nbr_list + cm_idx2;
   double lijsq1[num_nbr1];
   double lijsq2[num_nbr2];
   energy[0] = beobj.bending_energy_ipart(pos, nbrcm1,
                  num_nbr1, idx1, mesh.boxlen, btype1,
				//   mesh.lastbdry, mesh.pbc, 
				  lijsq1)
               +beobj.bending_energy_ipart(pos, nbrcm2,
                  	num_nbr2, idx2, mesh.boxlen, btype2,
					// mesh.lastbdry, mesh.pbc, 
					lijsq2);
   energy[0] += beobj.bending_energy_ipart_neighbour(pos, mesh, idx1)
               +beobj.bending_energy_ipart_neighbour(pos, mesh, idx2);
   energy[2] = chargeobj.debye_huckel_ipart(pos, idx1, mesh.N)
               +chargeobj.debye_huckel_ipart(pos, idx2, mesh.N);
   return energy[0]+energy[2];
}
// 
inline double McP::energy_mc_bechli(vector<double> &energy, Vec3d *pos, MESH_p mesh,
         	int idx1, int idx2, int cm_idx1, int cm_idx2, 
            int num_nbr1, int num_nbr2, BoundaryType btype1, BoundaryType btype2){
   double Etot = energy_mc_bech(energy, pos, mesh, idx1, idx2, cm_idx1,
               	cm_idx2, num_nbr1, num_nbr2, 
				      btype1, btype2);
   energy[4] = lineobj.energy_ipart(idx1) + lineobj.energy_ipart(idx2);
   return Etot+energy[4];
}
//
inline double McP::energy_mc_ch(vector<double> &energy, Vec3d *pos, MESH_p mesh,
      int idx1, int idx2, int cm_idx1, int cm_idx2, int num_nbr1, int num_nbr2){
   energy[2] = chargeobj.debye_huckel_ipart(pos, idx1, mesh.N)
               +chargeobj.debye_huckel_ipart(pos, idx2, mesh.N);
   return energy[2];
}
//
// int get_nstart(int lastbdry, int bdry_type){
//    int nframe;
//    if (bdry_type == 0 || bdry_type == 1) nframe = lastbdry+1;
//    else nframe = 0;
//    return nframe;
// }
//
int McP::monte_carlo_3d(Vec3d *pos, MESH_p mesh){
   int i, num_nbr, cm_idx;
   double x_o, y_o, z_o, x_n, y_n, z_n;
   double de, debe, dest, decharge, ders;
   double Einitot, Efintot;
   vector<double> Eini(5,0), Efin(5,0);
   double dxinc, dyinc, dzinc;
   double vol_i, vol_f;
   double dvol, de_vol, de_pressure;
   bool yes;
   int nframe;
   // the code does not generate random number for the boundary
   // if bdry_type == 0,1
   // FIX THIS BUG FOR FIXED BOUNDARY CONDITION
   // nframe = mesh.pbc ? 0 : (mesh.lastbdry + 1); 
   // nframe = get_nstart(mesh.lastbdry, steobj.bdry_type());
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   // std::cout << "RANK " << rank << " Algo pointer = " << (void*)(&Algo) << std::endl;
   nframe = 0;
   acceptedmoves = 0;
   if (!Algo) {
      std::cerr << "FATAL: Algo is null at start of monte_carlo_3d" << std::endl;
      abort();
   }
   for (i = 0; i < one_mc_iter; i++) {
      int idx = RandomGenerator::intUniform(nframe, mesh.N-1);
      cm_idx = idx*mesh.nghst;
      num_nbr = mesh.numnbr[idx];
      BoundaryType bt_idx = mesh.btype[idx];
      // cout << idx << " " << static_cast<int>(bt_idx) << endl;      
      Einitot = energy_mc_3d(Eini, pos, mesh, idx, mesh.nghst * idx, num_nbr, bt_idx);
      if (mesh.sphere){
         vol_i = steobj.volume_ipart(pos, (int *) (mesh.node_nbr_list + cm_idx),
               num_nbr, idx, mesh.boxlen, bt_idx
			);
      }
      x_o = pos[idx].x; y_o = pos[idx].y; z_o = pos[idx].z;
      //
      dxinc = (dfac) * (RandomGenerator::generateUniform(-1.0,1.0));
      dyinc = (dfac) * (RandomGenerator::generateUniform(-1.0,1.0));
      dzinc = (dfac) * (RandomGenerator::generateUniform(-1.0,1.0));
      //
      x_n = x_o + dxinc; y_n = y_o + dyinc; z_n = z_o + dzinc;
      pos[idx].x = x_n; pos[idx].y = y_n; pos[idx].z = z_n;
      Efintot = energy_mc_3d(Efin, pos, mesh, idx, mesh.nghst * idx, num_nbr, bt_idx);
		//
      de = Efintot - Einitot;
      if (mesh.sphere){
         vol_f = steobj.volume_ipart(pos,
               (int *) (mesh.node_nbr_list + cm_idx), num_nbr, idx,
               mesh.boxlen, bt_idx
			//    mesh.lastbdry, mesh.pbc
			);
         dvol = vol_f - vol_i;
         if(steobj.dovol()){
            de_vol = steobj.vol_energy_change(VolMonitored, dvol);
            de = de + de_vol;
         }
         if(steobj.dopressure()){
            de_pressure = steobj.PV_change(dvol);
            de = de + de_pressure;
         }
      }
      //
      if (!Algo) {
         std::cerr << "ERROR: Algo not initialized before MC step on rank "
                  << std::endl;
         abort();
      }
      if (!Algo) {
         std::cerr << "RANK " << rank << " Algo is EMPTY at MC step" << std::endl;
         abort();
      }
      yes = Algo(de, 0.0);
      if(yes){
         acceptedmoves +=  1;
         EneMonitored += de;
         bende += Efin[0]-Eini[0];
         stretche += Efin[1]-Eini[1];
         electroe += Efin[2]-Eini[2];
         linee += Efin[4]-Eini[4];
         vole += de_vol;
         VolMonitored += dvol;
      } else {
         pos[idx].x = x_o;
         pos[idx].y = y_o;
         pos[idx].z = z_o;
      }
  }
  return acceptedmoves;
}
//
int McP::monte_carlo_fluid(Vec3d *pos, MESH_p mesh){
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
  Vec3d aft_ij;

  double KAPPA;
  double av_bond_len=mesh.av_bond_len;
  bool yes, logic;

   // nframe = mesh.pbc ? 0 : (mesh.lastbdry + 1); if (nframe < 0) nframe = 0;
   // I may be wrong but I believe that boundary can be fluid too. 
   // even for pbc, the boundary can be fluid.
   nframe = 0;
   move = 0;

  int idxn, up, down;

  for (i = 0; i < one_mc_iter; i++) {
    // identify the pair to be divorced
    // stored as idx_del1 and idx_del2
    logic = false;
    while (!logic){
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
      }else {
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
//
inline double McP::energy_mc_be(vector<double> &energy, Vec3d *pos, 
        	MESH_p mesh, int idx1, int idx2, int cm_idx1, int cm_idx2, 
         	int num_nbr1, int num_nbr2,
			BoundaryType btype1, BoundaryType btype2){
   int *nbrcm1=mesh.node_nbr_list + cm_idx1;
   int *nbrcm2=mesh.node_nbr_list + cm_idx2;
   double lijsq1[num_nbr1];
   double lijsq2[num_nbr2];
   energy[0] = beobj.bending_energy_ipart(pos, nbrcm1,
                  num_nbr1, idx1, mesh.boxlen, btype1,
				//   mesh.lastbdry, mesh.pbc, 
				  lijsq1)
               +beobj.bending_energy_ipart(pos, nbrcm2,
                  num_nbr2, idx2, mesh.boxlen, btype2,
				//   mesh.lastbdry, mesh.pbc, 
				  lijsq2);
   energy[0] += beobj.bending_energy_ipart_neighbour(pos, mesh, idx1)
               +beobj.bending_energy_ipart_neighbour(pos, mesh, idx2);
   return energy[0];
}
//
int McP::local_idx(int num_nbr, int *node_nbr_list, int cm_idx){
   int idxn = RandomGenerator::intUniform(0, num_nbr-1);
   return node_nbr_list[cm_idx+idxn];
}
//
int McP::global_idx(int nframe, int N){
   return RandomGenerator::intUniform(nframe, N-1);
}
//
int McP::monte_carlo_lipid(Vec3d *pos, MESH_p mesh){
   int exchngdmoves = 0;
   int idx1, idx2, cm_idx1, cm_idx2;
   vector<double> Eini(5,0), Efin(5,0);
   bool yes, logic;
   int lip_idx1, lip_idx2, idxn, logic_break;
   double Einitot, Efintot, de;
   int num_nbr1, num_nbr2;
   int nframe=0;
   BoundaryType btype1;
   BoundaryType btype2;
   for (int i = 0; i < nexch_iter * one_mc_iter; ++i) {
      logic = true;
      idx1 = RandomGenerator::intUniform(nframe, mesh.N-1);
      cm_idx1 = mesh.nghst * idx1;
      num_nbr1 = mesh.numnbr[idx1];
      idx2 = get_idx2(num_nbr1, mesh.node_nbr_list, cm_idx1, nframe, mesh.N);
      cm_idx2 = mesh.nghst * idx2;
      num_nbr2 = mesh.numnbr[idx2];
      btype1 = mesh.btype[idx1];
      btype2 = mesh.btype[idx2];
      logic = (mesh.compA[idx1] == mesh.compA[idx2]);
      if (!logic){
         lip_idx1 = mesh.compA[idx1];
         lip_idx2 = mesh.compA[idx2];
         Einitot = energy_mc_exch(Eini, pos, mesh,
                                 idx1, idx2, 
                                 cm_idx1, cm_idx2,
                                 num_nbr1, num_nbr2,
                                 btype1, btype2);
         mesh.compA[idx2] = lip_idx1;
         mesh.compA[idx1] = lip_idx2;
         beobj.exchange(idx1, idx2, mesh);
         chargeobj.exchange(idx1,idx2);
         Efintot = energy_mc_exch(Efin, pos, mesh,
                                 idx1, idx2,
                                 cm_idx1, cm_idx2,
                                 num_nbr1, num_nbr2,
                                 btype1, btype2);
         de = Efintot-Einitot;
         yes = Boltzman(de, 0.0);
         if (yes){
            ++exchngdmoves;
            EneMonitored += de;
            bende += Efin[0]-Eini[0];
            stretche += Efin[1]-Eini[1];
            electroe += Efin[2]-Eini[2];
            linee += Efin[4]-Eini[4];
         }else{
            mesh.compA[idx1] = lip_idx1;
            mesh.compA[idx2] = lip_idx2;
            beobj.exchange(idx1, idx2, mesh);
            chargeobj.exchange(idx1,idx2);
         }
      }
   }
   return exchngdmoves;
}