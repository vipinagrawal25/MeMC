#ifndef CELLLIST_HPP
#define CELLLIST_HPP
#include <vector>
#include <array>
#include <cmath>
#include "mesh.hpp"
#include "vector.hpp"
// A class to handle the 3D cell list algorithm for molecular dynamics
class CellList {
public:
  void buildCellList(MESH_p mesh) {
    // Clear previous cell lists
    for (auto& cell : cellList_) {
        cell.clear();
    }

    // Assign particles to cells
    for (int i = 0; i < mesh.N; ++i) {
      
      int cellIdx = getCellIndex(mesh.pos[i].x, mesh.pos[i].y, mesh.pos[i].z, minL);
        // Remove the neighbours of the i here.
        cellList_[cellIdx].push_back(i);
    }
  }
  //
  void printParticlesInCells(Vec3d *particles_) {
    for (size_t cellIdx = 0; cellIdx < cellList_.size(); ++cellIdx) {
        std::cout << "Cell " << cellIdx << ": ";
        for (int pIdx : cellList_[cellIdx]) {
            std::cout << "(" << particles_[pIdx].x << ", " << particles_[pIdx].y << ", " << particles_[pIdx].z << ") ";
        }
        std::cout << std::endl;
    }
  }
  //
protected:
  double boxSize_;    // Size of the simulation box
  double cutoff_;     // Cutoff radius for interactions
  int numCells_;      // Number of cells per dimension
  double cellSize_;   // Size of each cell
  double minL;
  std::vector<std::vector<int>> cellList_;  // Cell list storing particle indices
  // std::vector<Particle> particles_;         // List of particles
  // Helper function to get the cell index based on particle position
  int getCellIndex(double x, double y, double z, double minL) {
    // Here need to add condition for x - minL - boxsize_. 
    // No particles can lie in the boundary
    // Else one edge of particles will repel the edge of other. 
    int ix = static_cast<int>((x - minL)/ cellSize_);
    int iy = static_cast<int>((y - minL)/ cellSize_);
    int iz = static_cast<int>((z - minL)/ cellSize_);
    // std::cout << x - minL << " " << y - minL << " " << cellSize_ << std::endl;
    return ix + numCells_ * (iy + numCells_ * iz);
  }
  // Get neighboring cells' indices (including itself)
  std::vector<int> getNeighboringCells(int cellIdx) {
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
  //
  bool isInArray(int arr[], int size, int value) {
    for (int i = 0; i < size; ++i) {
        if (arr[i] == value) {
            return true;
        }
    }
    return false;
  }
private:
};

#endif // CELLLIST_HPP