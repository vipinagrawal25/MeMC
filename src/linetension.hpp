#ifndef LINETENSION_HPP
#define LINETENSION_HPP
#include <string>
#include <vector>
#include <fstream>
#include "mesh.hpp"
#include "vector.hpp"
using namespace std;
class LTN {
public :
    LTN(const MESH_p& mesh, string fname);
    double energy_ipart(double *, int);
    double energy_ipart(int);
    double energy_total();
    bool calculate(){return lambda>0;}
private:
    double lambda;
    const MESH_p& mesh;
};
#endif