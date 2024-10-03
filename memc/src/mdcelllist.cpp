#include "mdcelllist.hpp"
#include <cmath>
#include <cstring>

// Constructor

extern "C" void  Cell_listread(bool *, double *, double *,
                               double *, double*, double *, char *);

MDCellList::MDCellList( string fname ) {
  char tmp_fname[128], temp_algo[128];
  string parafile, outfile;

  parafile = fname+"/para_file.in";
  sprintf(tmp_fname, "%s", parafile.c_str());
  Cell_listread(&doselfrepulsion, &boxSize_, &cutoff_,
                &minL, &sig, &epsl, tmp_fname);
  numCells_ = static_cast<int>(boxSize_ / cutoff_);
  cellSize_ = boxSize_ / numCells_;
  cellList_.resize(numCells_ * numCells_ * numCells_);
}

bool MDCellList::isSelfRepulsive(){
  return doselfrepulsion;
};

// int MDCellList::initCelllist(string fname) {
//   char tmp_fname[128], temp_algo[128];
//   string parafile, outfile;

//   parafile = fname+"/para_file.in";
//   sprintf(tmp_fname, "%s", parafile.c_str());
//   Cell_listread(&boxSize_, &cutoff_, &minL, tmp_fname);
//   numCells_ = static_cast<int>(boxSize_ / cutoff_);
//   // std::cout << boxSize_ << " " << cutoff_ << " " << numCells_ << " " << std::endl;
//   cellSize_ = boxSize_ / numCells_;
//   // numCells_ = 2;
//    cellList_.resize(numCells_ * numCells_ * numCells_);
//   // exit(0);
// }

// Helper function to get the 1D cell index for 3D space
int MDCellList::getCellIndex(double x, double y, double z, double minL) {
  // Here need to add condition for x - minL - boxsize_. No particles can lie in the boundary
  // Else one edge of particles will repel the edge of other. 
  int ix = static_cast<int>((x - minL)/ cellSize_);
  int iy = static_cast<int>((y - minL)/ cellSize_);
  int iz = static_cast<int>((z - minL)/ cellSize_);
  // std::cout << x - minL << " " << y - minL << " " << cellSize_ << std::endl;
    return ix + numCells_ * (iy + numCells_ * iz);
}

// Build the cell list by assigning particles to cells
void MDCellList::buildCellList(Vec3d *particles_, int N) {
    // Clear previous cell lists
    for (auto& cell : cellList_) {
        cell.clear();
    }

    // Assign particles to cells
    for (int i = 0; i < N; ++i) {
      
      int cellIdx = getCellIndex(particles_[i].x, particles_[i].y, particles_[i].z, minL);
        // Remove the neighbours of the i here.
        cellList_[cellIdx].push_back(i);
    }
}

// Get neighboring cells (including the cell itself)
std::vector<int> MDCellList::getNeighboringCells(int cellIdx) {
    std::vector<int> neighbors;

    int z = cellIdx / (numCells_ * numCells_);
    int y = (cellIdx % (numCells_ * numCells_)) / numCells_;
    int x = cellIdx % numCells_;

    // Loop through neighboring cells in 3x3x3 cube
    for (int dz = -1; dz <= 1; ++dz) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                int nx = (x + dx + numCells_) % numCells_;
                int ny = (y + dy + numCells_) % numCells_;
                int nz = (z + dz + numCells_) % numCells_;
                neighbors.push_back(nx + numCells_ * (ny + numCells_ * nz));
            }
        }
    }

    return neighbors;
}

bool isInArray(int arr[], int size, int value) {
    for (int i = 0; i < size; ++i) {
        if (arr[i] == value) {
            return true;
        }
    }
    return false;
}

double MDCellList::totalRepulsiveEnergy(Vec3d *particles_, MESH_p mesh){
  int i; 
  double et = 0.0;
  for(i=0; i< mesh.N; i++){
    et += computeSelfRep(particles_, mesh, i);
  }
  return et;
}

// Compute forces between particles using the cell list
double MDCellList::computeSelfRep(Vec3d *particles_, MESH_p mesh, int idx) {

  int cellIdx = getCellIndex(particles_[idx].x, particles_[idx].y, particles_[idx].z, minL);
  int num_nbr = mesh.numnbr[idx];
  int cm_idx = idx*mesh.nghst;

  // for (size_t cellIdx = 0; cellIdx < cellList_.size(); ++cellIdx) {
  std::vector<int> neighbors = getNeighboringCells(cellIdx);
  Vec3d p1 = particles_[idx];
  double EselfRep = 0.0;
  double sigr2 = 0.0;

  // std::cout << particles_[idx].x << " " << particles_[idx].y << " " << particles_[idx].z << std::endl;
  for (int j : cellList_[cellIdx]) {
    if(isInArray((int *)(mesh.node_nbr_list + cm_idx), num_nbr, j)){
      // std::cout << particles_[i].x << " " << particles_[i].y << " " << particles_[i].z << std::endl;
        Vec3d p2 = particles_[j];
        double dx = p1.x - p2.x;
        double dy = p1.y - p2.y;
        double dz = p1.z - p2.z;
        double r2 = dx * dx + dy * dy + dz * dz;
        sigr2 = sig*sig/r2;
        EselfRep += (4*epsl*std::pow(sigr2, 3)); // Example Lennard-Jones force

    } 
    // std::cout << particles_[i].x << " " << particles_[i].y << " " << particles_[i].z << std::endl;
  }
  // loop over all the neighbors
  for (int neighborIdx : neighbors) {
    for (int j : cellList_[neighborIdx]) {
      if(isInArray((int *)(mesh.node_nbr_list + cm_idx), num_nbr, j)){
        Vec3d p2 = particles_[j];
        // std::cout << particles_[j].x << " " << particles_[j].y << " " << particles_[j].z << std::endl;
        // Calculate distance between particles
        double dx = p1.x - p2.x;
        double dy = p1.y - p2.y;
        double dz = p1.z - p2.z;
        double r2 = dx * dx + dy * dy + dz * dz;
        sigr2 = sig*sig/r2;
        EselfRep += (4*epsl*std::pow(sigr2, 3)); // Example Lennard-Jones force
      }
    }
  }
  return EselfRep;
}
    

// }

// Print the particles in each cell for debugging
void MDCellList::printParticlesInCells(Vec3d *particles_) {
    for (size_t cellIdx = 0; cellIdx < cellList_.size(); ++cellIdx) {
        std::cout << "Cell " << cellIdx << ": ";
        for (int pIdx : cellList_[cellIdx]) {
            std::cout << "(" << particles_[pIdx].x << ", " << particles_[pIdx].y << ", " << particles_[pIdx].z << ") ";
        }
        std::cout << std::endl;
    }
}

