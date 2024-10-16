#include "selfavoidance.hpp"

extern "C" void  SelfAvoidRead(bool *, double *, double *,
                               double *, double*, double *, char *);

SelfAvoid::SelfAvoid(MESH_p mesh, string fname){
  char tmp_fname[128];
  string parafile, outfile;

  parafile = fname+"/para_file.in";
  sprintf(tmp_fname, "%s", parafile.c_str());
  SelfAvoidRead(&doselfrepulsion, &boxSize_, &cutoff_,
                &minL, &sig, &epsl, tmp_fname);
  numCells_ = static_cast<int>(boxSize_ / cutoff_);
  cellSize_ = boxSize_/numCells_;
  cellList_.resize(numCells_ * numCells_ * numCells_);
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
  double sigr2 = 0.0;

  // std::cout << particles_[idx].x << " " << particles_[idx].y << " " << particles_[idx].z << std::endl;
  for (int j : cellList_[cellIdx]){
    if(isInArray((int *)(mesh.node_nbr_list + cm_idx), num_nbr, j)){
      // std::cout << particles_[i].x << " " << particles_[i].y << " " << particles_[i].z << std::endl;
        Vec3d p2 = mesh.pos[j];
        double dx = p1.x - p2.x;
        double dy = p1.y - p2.y;
        double dz = p1.z - p2.z;
        double r2 = dx * dx + dy * dy + dz * dz;
        sigr2 = sig*sig/r2;
        EselfRep += (4*epsl*std::pow(sigr2, 6));
        // Example Lennard-Jones force
    } 
  }
  // loop over all the neighbors
  for (int neighborIdx : neighbors) {
    for (int j : cellList_[neighborIdx]) {
      if(isInArray((int *)(mesh.node_nbr_list + cm_idx), num_nbr, j)){
        Vec3d p2 = mesh.pos[j];
        double dx = p1.x - p2.x;
        double dy = p1.y - p2.y;
        double dz = p1.z - p2.z;
        double r2 = dx * dx + dy * dy + dz * dz;
        sigr2 = sig*sig/r2;
        EselfRep += (4*epsl*std::pow(sigr2, 6)); // Example Lennard-Jones force
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
  return et;
}