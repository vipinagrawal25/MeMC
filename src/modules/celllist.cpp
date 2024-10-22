#include "celllist.hpp"
#include <cmath>
#include <cstring>

// Constructor

extern "C" void  Cell_listread(bool *, double *, double *,
                               double *, double*, double *, char *);

// CellList::CellList( string fname ) {
//   char tmp_fname[128], temp_algo[128];
//   string parafile, outfile;

//   parafile = fname+"/para_file.in";
//   sprintf(tmp_fname, "%s", parafile.c_str());
//   Cell_listread(&doselfrepulsion, &boxSize_, &cutoff_,
//                 &minL, &sig, &epsl, tmp_fname);
//   numCells_ = static_cast<int>(boxSize_ / cutoff_);
//   cellSize_ = boxSize_ / numCells_;
//   cellList_.resize(numCells_ * numCells_ * numCells_);
// }

// bool CellList::isSelfRepulsive(){
//   return doselfrepulsion;
// };

// int CellList::initCelllist(string fname) {
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

