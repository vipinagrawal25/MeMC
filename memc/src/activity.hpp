#ifndef ACTIVITY_HPP
#define ACTIVITY_HPP

#include <string>
#include <vector>
#include "random_gen.hpp"
// #include "mesh.hpp"
#include "vector.hpp"

class ACT {
public : 
  int initACT(int , std::string );
  double getActivityIdx(int );
  bool is_active();
private:
    bool doactivity;  
    vector <double> Actarr;
};


#endif
