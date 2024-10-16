#ifndef BENDING_HPP
#define BENDING_HPP
#include <string>
#include "mesh.hpp"
#include "vector.hpp"
#include "misc.hpp"
#include <vector>
// #include "global.h"
class BE{
public :
  BE(const MESH_p& mesh, std::string fname);
  double bending_energy_ipart_neighbour(Vec3d *pos, MESH_p mesh, int idx);
  double bending_energy_ipart(Vec3d *pos, int *node_nbr, int num_nbr,
                              int idx, int, double, int);
  double bending_energy_total(Vec3d *pos, MESH_p mesh);
  void init_coefbend(int *lipA, int N);
  void printbend(){print(coef_bend);}
  int getbend(int i){return coef_bend[i];}
  double getbend1(){return bend1;}
  double getbend2(){return bend2;}
  void exchange(int idx1, int idx2);
private:
  double bend1, bend2, spC1, spC2;
  std::vector<double> coef_bend;
  // std::vector<double> spcurv;
  double spcurv;
/*------------------------*/
double cotangent(Vec3d si, Vec3d sk, Vec3d sj){
  ///
  ///  @param si  coordinate of ith point
  ///  @param sk  coordinate of kth point
  ///  @param sj  coordinate of jth point 
   
  ///  @return   ({si-sk}.{sj-sk})/sqrt(({si-sk}x{sj-sk})^2)
  /// angle between vector si and sj
  ///
  Vec3d drik, drjk, cross;
  double cot_theta;  
  double inner_prod;
  //
  drik = si - sk; 
  drjk = sj - sk; 
  cross = cross_product(drik, drjk);
  inner_prod = inner_product(drik, drjk);
  cot_theta = inner_prod/sqrt(inner_product(cross,cross));
  //
  return cot_theta;
}
/*------------------------*/
double cotangent(double a, double b, double c){
   /// @brief  a, b, c are the length of the sides of a triangle 
   ///  @return   0.25*(a*a+b*b-c*c)/area; where area is the area of triangle
  double s = 0.5*(a+b+c);
  double area = sqrt(s*(s-a)*(s-b)*(s-c));
  double cot_theta=0.25*(a*a+b*b-c*c)/area;
  return cot_theta;
}
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