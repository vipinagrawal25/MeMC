#ifndef FILE_Position_SEEN
#define FILE_Position_SEEN
/*-----------------------------------------*/
#include <iostream>
#include <math.h>
using namespace std;
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
class Vec2d{
public:
  double x,y;
  Vec2d(double a, double b){
     x=a, y=b;
  }
  Vec2d();
private:
};

inline Vec2d::Vec2d(){
  x = 0.0;
  y = 0.0;  
}
/*-----------------------------------------*/
inline Vec3d::Vec3d(){
  x = 0.0;
  y = 0.0;
  z = 0.0;
}
/*---------------------------------------*/
inline Vec3d::Vec3d(int a, int b, int c){
  x = (double) a;
  y = (double) b;
  z = (double) c;
}
/*---------------------------------------*/
inline Vec3d::Vec3d(float a, float b, float c){
  x = (double) a;
  y = (double) b;
  z = (double) c;
}
/*---------------------------------------*/
inline Vec3d::Vec3d(double a, double b, double c){
  x =  a;
  y =  b;
  z =  c;
}
/*---------------------------------------*/
inline Vec3d Vec3d::operator+(Vec3d param){
  Vec3d temp;
  temp.x = x+param.x;
  temp.y = y+param.y;
  temp.z = z+param.z;
  return(temp);
}
/*---------------------------------------*/
inline Vec3d Vec3d::operator-(Vec3d param){
  Vec3d temp;
  temp.x = x-param.x;
  temp.y = y-param.y;
  temp.z = z-param.z;
  return(temp);
}
/*---------------------------------------*/
inline Vec3d Vec3d::operator*(double param){
  Vec3d temp;
  temp.x=param*x;
  temp.y=param*y;
  temp.z=param*z;
  return(temp);
}
/*---------------------------------------*/
inline void print(Vec3d a){
  cout<<a.x<<"\t"<<a.y<<"\t"<<a.z<<"\n";
}
/*---------------------------------------*/
inline void print(Vec2d a){
  cout<<a.x<<"\t"<<a.y<<"\n";
}
/*---------------------------------------*/
inline void print(Vec3d* a, int N){
  for (int i = 0; i < N; ++i){
    cout<<a[i].x<<"\t"<<a[i].y<<"\t"<<a[i].z<<"\n";
  }
}
/*---------------------------------------*/
inline Vec3d Vec3d::operator/(double param){
  Vec3d temp;
  temp.x=x/param;
  temp.y=y/param;
  temp.z=z/param;
  return(temp);
}
/*---------------------------------------*/
inline void printnorm(Vec3d* a, int N){
  for (int i = 0; i < N; ++i){
    cout<<sqrt(a[i].x*a[i].x + a[i].y*a[i].y + a[i].z*a[i].z)<<"\n";
  }
}
/*---------------------------------------*/
inline double inner_product(Vec3d s1, Vec3d s2){
    return s1.x*s2.x + s1.y*s2.y + s1.z*s2.z;
}
/*---------------------------------------*/
inline double norm(Vec3d s1){
    return sqrt(inner_product(s1,s1));
}
/*---------------------------------------*/
inline double normsq(Vec3d s1){
    return inner_product(s1,s1);
}
/*---------------------------------------*/
inline Vec3d cross_product(Vec3d s1, Vec3d s2){
  Vec3d crosprod;
  crosprod.x = s1.y*s2.z - s1.z*s2.y;
  crosprod.y = s1.z*s2.x - s1.x*s2.z;
  crosprod.z = s1.x*s2.y - s1.y*s2.x;
  return crosprod;
}
/*---------------------------------------*/
#endif /* !FILE_Position_SEEN */
