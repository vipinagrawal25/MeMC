#ifndef BENDING_HPP
#define BENDING_HPP
#include <string>
#include "mesh.hpp"
#include "vector.hpp"
#include "misc.hpp"
#include <vector>
#include <functional>
// #include "global.h"
class BE{
public :
  BE(const MESH_p& mesh, std::string fname);
  double bending_energy_ipart_neighbour(Vec3d *pos, MESH_p mesh, int idx);
  // double bending_energy_ipart(Vec3d *pos, int *node_nbr, int num_nbr,
  //     int, int, double, int);
  // double bending_energy_ipart(Vec3d *pos, int *node_nbr, int num_nbr,
  //     int, int, double, int, double*);
  double bending_energy_total(Vec3d *pos, MESH_p mesh);
  void printbend(const MESH_p& mesh){
    for (int i = 0; i < mesh.nghst * mesh.N; ++i){
      cout << bendij[i] << " ";
    }
  }
  int getbend(int i){return coef_bend[i];}
  double getbend1(){return bend1;}
  double getbend2(){return bend2;}
  std::function<double(Vec3d *, int *, int , int, int , double , int ,
          double *)> bending_energy_ipart;
  double SeungNelson(Vec3d *, int *, int , int, int , double , int , double *);
  double Itzykson(Vec3d *, int *, int , int, int , double , int , double *);
  bool isexch(){return multicomp;}
  std::function<void(int, int, const MESH_p&)> exchange;
  void set_nbrbending(int idx1, const MESH_p& mesh);
private:
  double bend1, bend2, spC1, spC2;
  std::vector<double> coef_bend;
  std::vector<double> bendij;
  double spcurv;
  bool iGauss;
  bool multicomp;
  string method="SN";
  // This is required by SeungNelson. Right way is to pass it, but I do not want
  // to change the function call everywhere.
  int ghost;
  void init_coefbend(int *lipA, int N);
  void init_bendij(MESH_p mesh);
  void exchange_node(int idx1, int idx2);
  void exchange_bond(int idx1, int idx2, const MESH_p& mesh);
/*------------------------*/
double voronoi_area(double cotJ, double cotK, 
      double jsq, double ksq, double area){
  /// @brief Estimate the area of the voronoi cell. If I J K are the nodes of
  /// triangle
  ///  @param cotJ angle at node j (see paper/paper.pdf)
  ///  @param cotK angle at node k (see paper/paper.pdf)
  ///  @param jsqr square of the length of bond i-k
  /// @param ksq square of the length of bond i-j  
  ///  @param area area of the triangle 
  /// @return  Given two cotangent angles, it returns either the area due to perpendicular bisector,
  /// or the barycenter.
 double sigma;
  if (cotJ>0 && cotK>0){
      if (cotJ*cotK<1){
          // all angles are acute;
          sigma = 0.125*(cotJ*jsq+cotK*ksq);
      }else{
          sigma = 0.5*area;
      }
  }else{
     sigma = 0.25*area;
 }
  return sigma;
}
/*-------------------------------------------------*/
Vec3d diff(Vec3d a, Vec3d b, double lenth){
   Vec3d ab=a-b;
   if (ab.x >= 0.5 * lenth) ab.x -= lenth;
   if (ab.x < -0.5 * lenth) ab.x += lenth;

   if (ab.y >= 0.5 * lenth) ab.y -= lenth;
   if (ab.y < -0.5 * lenth) ab.y += lenth;

   if (ab.z >= 0.5 * lenth) ab.z -= lenth;
   if (ab.z < -0.5 * lenth) ab.z += lenth;

   return ab;
}
/*-------------------------------------------------*/
};
#endif