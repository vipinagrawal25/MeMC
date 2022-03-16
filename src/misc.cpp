#include <iostream>
#include <sys/stat.h>
#include "math.h"
#include <string>
#include "../include/Position.h"
using namespace std;
/*-----------------------------------------------*/
double SqEr(double Arr1[], double Arr2[],int nn){
  double error=0;
  for (int cdim=0; cdim < nn; cdim++){
    error=error+(Arr2[cdim]-Arr1[cdim])*(Arr2[cdim]-Arr1[cdim]);
  }
  error = sqrt(error);
  return error;
}
/*-----------------------------------------------*/
bool FileExists(const std::string &s){
  struct stat buffer;
  return (stat (s.c_str(), &buffer) == 0);
}
/*-----------------------------------------------*/
template<typename T>
T absolute(T value){
  if (value>=0){
    return value;
  }else{
    return value*(-1);
  }
}
/*----------------------------------------------*/
void __attribute__((weak)) check_param(){
  cout << "I believe all your model parameters are physical. Otherwise, define function: "
          "void check_param() in model.cpp file" << endl;
}
// -----------------------------------------------
void __attribute__((weak)) write_param(string fname){
  cout << "I can not find any implementation to write model parameters." 
          "Hence I will not write anything." << endl;
}
/*-----------------------------------------------*/
void print(double *arr, int start, int skip, int end){
  for (int in = 0; in < int(end/skip); ++in){
      cout << arr[skip*in+start] << " ";
  }
  cout << endl;
}
/*-----------------------------------------------*/
void print(double *arr, int nn){
  for (int in = 0; in < nn; ++in){
      cout << arr[in] << "\t";
  }
  cout << endl;
}
/*-----------------------------------------------*/
void print(double *arr, int nn, int pp){
  for (int in = 0; in < nn; ++in){
    for (int ip = 0; ip < pp; ++ip){
      cout << arr[in*pp+ip] << "\t";
    }
    cout << endl;
  }
}
/*-----------------------------------------------*/
// void print(vec2 *arr, int nn){
//   for (int in = 0; in < nn; ++in){
//       cout << arr[in].x << "\t";
//       cout << arr[in].y << "\n";
//   }
// }
/*-----------------------------------------------*/
void add(double *added, double *arr1, double *arr2, int nn){
  for (int ii = 0; ii < nn; ++ii){
    added[ii] = arr1[ii] + arr2[ii];
  }
}
/*-----------------------------------------------*/
void substract(double *subsed, double *arr1, double *arr2, int nn ){
  for (int ii = 0; ii < nn; ++ii){
    subsed[ii] = arr1[ii] - arr2[ii];
  }
}
/*-----------------------------------------------*/
double norm(double *y, int nn ){
  double norr =0;
  for (int ii = 0; ii < nn; ++ii){
    norr = norr+y[ii]*y[ii];
  }
  return sqrt(norr);
}
/*-----------------------------------------------*/
double norm(const double *y, int nn ){
  double norr =0;
  for (int ii = 0; ii < nn; ++ii){
    norr = norr+y[ii]*y[ii];
  }
  return sqrt(norr);
}
/*-----------------------------------------------*/
// template<typename T, typename ... Args>
// T First(T first,Args ... args ){
//   return first;
// }
// /*-----------------------------------------------*/
// template<typename T, typename ... Args>
// T Last(T first, T second, Args ... args ){
//   return Last(args...);
// }
// /*-----------------------------------------------*/
// double upScale(double *y, int factor, int nn ){
//   for (int in = 0; in < nn; ++in){
//     y[in] = y[in]*factor;
//   }
// }
/*-----------------------------------------------*/
void downScale(double *yscaled, double *y, int factor, int nn ){
  for (int in = 0; in < nn/factor; ++in){
    yscaled[in] = y[in*factor];
  }
}
/*-----------------------------------------------*/
void zeros(double *yzero, int ndim){
  for (int idim = 0; idim < ndim; ++idim){
    yzero[idim]=0.;
  }
}
/*-----------------------------------------------*/
inline double pos_coord(POSITION pos,char dirn='z'){
  if (dirn=='x'){return pos.x;}
  else if(dirn=='y'){return pos.y;}
  else if(dirn=='z'){return pos.z;}
}
/*-----------------------------------------------*/
void max(int *amaxind, double *amaxval, POSITION *pos, int ndim,char dirn){
  // function returns the value and index of the maximum entry.
  int maxind=0;
  double maxval=-1e+16;
  //
  for (int idim = 0; idim < ndim; ++idim){
    if(pos_coord(pos[idim],dirn)>maxval){
      maxval = pos_coord(pos[idim],dirn);
      maxind = idim;
    }
  }
  *amaxind=maxind;
  *amaxval=maxval;
}
/*-----------------------------------------------*/
void min(int *aminind, double *aminval, POSITION *pos, int ndim,char dirn){
  // function returns the value and index of the maximum entry.
  int minind=0;
  double minval=1e+16;
  //
  for (int idim = 0; idim < ndim; ++idim){
    if(pos_coord(pos[idim],dirn)<minval){
      minval = pos_coord(pos[idim],dirn);
      minind = idim;
    }
  }
  *aminind=minind;
  *aminval=minval;
}
/*-----------------------------------------------*/