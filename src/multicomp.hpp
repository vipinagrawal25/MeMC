#ifndef MULTICOMP_HPP
#define MULTICOMP_HPP
#include "misc.hpp"
#include "mesh.hpp"
#include "vector.hpp"
#include <string>
#include <vector>
#include "bending.hpp"
class MulCom{
public:
   MulCom(const MESH_p& mesh, std::string fname);
   Vec2d reg_soln_ipart(Vec3d *pos, MESH_p mesh, int idx);
   double reg_soln_ipart_neighbour(Vec3d *pos, MESH_p mesh, int idx);
   double reg_soln_tot(Vec3d *pos, MESH_p mesh);
   double getepssqby2(){return epssqby2;}
   double gradphisq_ipart(Vec3d *pos, MESH_p mesh, int idx);
   double gradphisq_ipart_neighbour(Vec3d *pos, MESH_p mesh, int idx);
   bool calculate(){return iregsoln;}
private:
   // MESH_p &mesh;
   bool iregsoln;
   double kai;
   double epssqby2;
   double phi_ipart(int *lipA, int *node_nbr, int num_nbr, int idx);
   void phi_ipart_neighbour(double *phi, MESH_p mesh, int idx);
   Vec2d gradphisq(double *phi, Vec3d *pos, int *node_nbr, int num_nbr, int idx,
        int bdry_type, double lenth, int edge);
};
#endif