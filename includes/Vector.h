#ifndef FILE_Position_SEEN
#define FILE_Position_SEEN
/*-----------------------------------------*/
#include <iostream>
#include <math.h>
using namespace std;
/*----------------------------------------*/
// typedef class Vec3d Vec3d;
/*-----------------------------------------*/
typedef struct{
  double x,y;
}Vec2d;
/*-----------------------------------------*/
class Vec3d{
public:
  double x,y,z;
  Vec3d();
  Vec3d(int,int,int);
  Vec3d(float,float,float);
  Vec3d(double,double,double);
  // definining operators 
  Vec3d operator+(Vec3d);
  Vec3d operator-(Vec3d);
  Vec3d operator*(double);
  Vec3d operator/(double);
  // Vec3d (double)*operator;
private:
};
/*-----------------------------------------*/
inline Vec3d::Vec3d(){
  x = 0.0;
  y = 0.0;
  z = 0.0;
}
inline Vec3d::Vec3d(int a, int b, int c){
  x = (double) a;
  y = (double) b;
  z = (double) c;
}
inline Vec3d::Vec3d(float a, float b, float c){
  x = (double) a;
  y = (double) b;
  z = (double) c;
}
inline Vec3d::Vec3d(double a, double b, double c){
  x =  a;
  y =  b;
  z =  c;
}
inline Vec3d Vec3d::operator+(Vec3d param){
  Vec3d temp;
  temp.x = x+param.x;
  temp.y = y+param.y;
  temp.z = z+param.z;
  return(temp);
}
inline Vec3d Vec3d::operator-(Vec3d param){
  Vec3d temp;
  temp.x = x-param.x;
  temp.y = y-param.y;
  temp.z = z-param.z;
  return(temp);
}
inline Vec3d Vec3d::operator*(double param){
  Vec3d temp;
  temp.x=param*x;
  temp.y=param*y;
  temp.z=param*z;
  return(temp);
}
inline void print(Vec3d a){
  cout<<a.x<<"\t"<<a.y<<"\t"<<a.z<<"\n";
}
/*---------------------------------------*/
inline Vec3d Vec3d::operator/(double param){
  Vec3d temp;
  temp.x=x/param;
  temp.y=y/param;
  temp.z=z/param;
  return(temp);
}
double inner_product(Vec3d s1, Vec3d s2);
double norm(Vec3d s1);
double normsq(Vec3d s1);
Vec3d Vec3d_add(Vec3d s1, Vec3d s2, double fac);
    /* Returns s1 + fac*s2*/
Vec3d cross_product(Vec3d s1, Vec3d s2);
Vec3d diff_pbc(Vec3d r1, Vec3d r2, double len);
#endif /* !FILE_Position_SEEN */