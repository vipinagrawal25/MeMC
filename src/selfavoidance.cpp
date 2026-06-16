#include "selfavoidance.hpp"

extern "C" void  SelfAvoidRead(bool *, double *, double *, double *, double*, char *);

SelfAvoid::SelfAvoid(MESH_p mesh, string fname){
  char tmp_fname[128];
  string parafile, outfile;
  double cutoff_scale, minL_scale;

  parafile = fname+"/para_file.in";
  sprintf(tmp_fname, "%s", parafile.c_str());
  SelfAvoidRead(&doselfrepulsion, &cutoff_scale, &minL_scale, &sig_scale, &epsl_scale, tmp_fname);
  
  // Store av_bond_len as scale reference
  av_bond_len_scale = mesh.av_bond_len;
  
  // Apply scale factors to convert from scale factors to actual values
  cutoff_ = cutoff_scale * av_bond_len_scale;
  // minL is the minimum coordinate - scale it as a fraction of the negative radius
  // minL_scale of 0.5 means minL = -0.5*mesh.radius
  minL = -mesh.radius * minL_scale;
  
  boxSize_ = 2*(mesh.radius-minL);
  numCells_ = static_cast<int>(boxSize_ / cutoff_);
  cellSize_ = boxSize_/numCells_;
  
  // Initialize kBT to 1.0 (will be set properly by setKBT)
  kBT = 1.0;
  // cellList_.resize(numCells_ * numCells_ * numCells_);
}
//
double SelfAvoid::LJ(Vec3d p1, Vec3d p2){
  // Vec3d p2 = mesh.pos[j];
  double dx = p1.x - p2.x;
  double dy = p1.y - p2.y;
  double dz = p1.z - p2.z;
  double r2 = dx * dx + dy * dy + dz * dz;
  double sig_actual = sig_scale * av_bond_len_scale;
  double sigr2 = sig_actual*sig_actual/r2;
  double epsl = epsl_scale * kBT;
  // 
  double en=4*epsl*std::pow(sigr2, 6);
  return en;
}
//
double SelfAvoid::computeSelfRep(MESH_p mesh, int idx) {
  int cellIdx = getCellIndex(mesh.pos[idx].x, mesh.pos[idx].y, mesh.pos[idx].z, 
                minL);
  int num_nbr = mesh.numnbr[idx];
  int cm_idx = idx*mesh.nghst;

  // for (size_t cellIdx = 0; cellIdx < cellList_.size(); ++cellIdx) {
  std::vector<int> neighbors = getNeighboringCells(cellIdx);
  Vec3d p1 = mesh.pos[idx];
  double EselfRep = 0.0;
  // double sigr2 = 0.0;

  // std::cout << particles_[idx].x << " " << particles_[idx].y << " "
  //  << particles_[idx].z << std::endl;
//   for (int j : cellList_[cellIdx]){
//     if(isInArray((int *)(mesh.node_nbr_list + cm_idx), num_nbr, j)){
//         EselfRep += LJ(p1, mesh.pos[j]);
//      } 
//   }
  // loop over all the neighbors
  for (int neighborIdx : neighbors) {
    for (int j : cellList_[neighborIdx]){
      if(j != idx && isNotInArray((int *)(mesh.node_nbr_list + cm_idx), num_nbr, j)){
        EselfRep += LJ(p1, mesh.pos[j]);
      }
    }
  }
  return EselfRep;
}

double SelfAvoid::totalRepulsiveEnergy(MESH_p mesh){
  int i; 
  double et = 0.0;
  for(i=0; i< mesh.N; i++){
    et += computeSelfRep(mesh, i);
  }
  return 0.5*et;
}