#ifndef FILE_Position_SEEN
#define FILE_Position_SEEN
/*-----------------------------------------*/
#include <iostream>
#include <math.h>
using namespace std;
/*----------------------------------------*/
// typedef class POSITION POSITION;
/*-----------------------------------------*/
class POSITION{
public:
  double x,y,z;
  POSITION();
  POSITION(int,int,int);
  POSITION(float,float,float);
  POSITION(double,double,double);
  // definining operators 
  POSITION operator+(POSITION);
  POSITION operator-(POSITION);
  POSITION operator*(double);
  POSITION operator/(double);
  // POSITION (double)*operator;
private:
};
/*-----------------------------------------*/
inline POSITION::POSITION(){
  x = 0.0;
  y = 0.0;
  z = 0.0;
}
inline POSITION::POSITION(int a, int b, int c){
  x = (double) a;
  y = (double) b;
  z = (double) c;
}
inline POSITION::POSITION(float a, float b, float c){
  x = (double) a;
  y = (double) b;
  z = (double) c;
}
inline POSITION::POSITION(double a, double b, double c){
  x =  a;
  y =  b;
  z =  c;
}
inline POSITION POSITION::operator+(POSITION param){
  POSITION temp;
  temp.x = x+param.x;
  temp.y = y+param.y;
  temp.z = z+param.z;
  return(temp);
}
inline POSITION POSITION::operator-(POSITION param){
  POSITION temp;
  temp.x = x-param.x;
  temp.y = y-param.y;
  temp.z = z-param.z;
  return(temp);
}
inline POSITION POSITION::operator*(double param){
  POSITION temp;
  temp.x=param*x;
  temp.y=param*y;
  temp.z=param*z;
  return(temp);
}
inline double inner_product(POSITION s1, POSITION s2){
    return s1.x*s2.x + s1.y*s2.y + s1.z*s2.z;
}
inline double norm(POSITION s1){
    return sqrt(inner_product(s1,s1));
}
inline double normsq(POSITION s1){
    return inner_product(s1,s1);
}
inline POSITION Position_add(POSITION s1, POSITION s2, double fac){
    /* Returns s1 + fac*s2*/
    POSITION add;
    add.x = s1.x + fac*s2.x;
    add.y = s1.y + fac*s2.y;
    add.z = s1.z + fac*s2.z;
    return add;
}
inline POSITION cross_product(POSITION s1, 
        POSITION s2){

    POSITION crosprod;
    crosprod.x = s1.y*s2.z - s1.z*s2.y;
    crosprod.y = s1.z*s2.x - s1.x*s2.z;
    crosprod.z = s1.x*s2.y - s1.y*s2.x;
    return crosprod;

}
inline POSITION POSITION::operator/(double param){
  POSITION temp;
  temp.x=x/param;
  temp.y=y/param;
  temp.z=z/param;
  return(temp);
}
/*---------------------------------------*/
inline void print(POSITION a){
  cout<<a.x<<"\t"<<a.y<<"\t"<<a.z<<"\n";
}
/*---------------------------------------*/
#endif /* !FILE_Position_SEEN */
