#include <iostream>
#include <sys/stat.h>
#include "math.h"
#include <string>
#include "Vector.h"
#include "global.h"
#include <iomanip>
#include <sstream>
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
void print(float *arr, int start, int skip, int end){
  for (int in = 0; in < int(end/skip); ++in){
      cout << arr[skip*in+start] << " ";
  }
  cout << endl;
}
/*-----------------------------------------------*/
void print(double *arr, int nn){
  for (int in = 0; in < nn; ++in){
      cout << arr[in] << "\n";
  }
  cout << endl;
}
/*-----------------------------------------------*/
void print(double *arr, int nn, string fname){
  ofstream outfile;
  outfile.open(fname);
  for (int in = 0; in < nn; ++in){
      outfile << arr[in] << "\n";
  }
  outfile.close();
}
/*-----------------------------------------------*/
void print(int *arr, int nn){
  for (int in = 0; in < nn; ++in){
      cout << arr[in] << "\n";
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

inline double pos_coord(Vec3d pos, char dirn='z'){
  if (dirn=='x'){return pos.x;}
  else if(dirn=='y'){return pos.y;}
  else if(dirn=='z'){return pos.z;}
  return 0e0;
}
/*-----------------------------------------------*/
void max(int *amaxind, double *amaxval, Vec3d *pos, int ndim,char dirn){
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
void min(int *aminind, double *aminval, Vec3d *pos, int ndim,char dirn){
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
/*-----------------------------------------------*/
double height_rms(Vec3d *Pos, MBRANE_p mbrane){
  double radius=mbrane.radius;
  double N=mbrane.N;
  double hrms=0;
  double hh;
  for (int i = 0; i < N; ++i){
    hh = sqrt(Pos[i].x*Pos[i].x+Pos[i].y*Pos[i].y+Pos[i].z*Pos[i].z) - radius;
    hrms += hh*hh;
  }
  return sqrt(hrms/N);
}
/*-----------------------------------------------*/
int get_nstart(int N, int bdrytype){
    static int nf1;
    int nf2;
    nf1 = (int) sqrt((double)N);
    switch (bdrytype) {
        case 0:
            nf2 = 2 * nf1;
            break;
        case 1:
            nf2 = 4 * nf1; 
            break;
        default:
            nf2 = 0;
    }
    return nf2;;
}
/*-----------------------------------------------*/
// just in case I need it later;
int print_sanity(int *nbr_del1, int *nbr_del2, int *nbr_add1, int *nbr_add2,
                 int del1, int del2, int add1, int add2, char *fname) {
  FILE *fid;
  int k = 10;

  fid = fopen(fname, "a");
  /* fprintf(fid, "start\n"); */
  fprintf(fid, "%06d\t%06d\t%06d\t%06d\t%06d\n", k, del1, del2, add1, add2);
  for (int j = 0; j < 12; j++) {
    fprintf(fid, "%06d\t%06d\t%06d\t%06d\t%06d\n", j, nbr_del1[j], nbr_del2[j],
            nbr_add1[j], nbr_add2[j]);
  }
  fprintf(fid, "next line\n");
  fflush(fid);
  fclose(fid);
  return 0;
}

int print_sanity(Vec3d *pos, int *nbr_del1, int *nbr_del2, int *nbr_add1,
                 int *nbr_add2, int del1, int del2, int add1, int add2,
                 char *ftag, int idx) {
  FILE *fid;
  char fname[128];
  int k = 10;

  sprintf(fname, "check_pos/%s_pos_%04d_.dat", ftag, idx);
  fid = fopen(fname, "w");
  /* fprintf(fid, "start\n"); */
  fprintf(fid, "1 %g\t%g\n", pos[del1].x, pos[del1].y);
  for (int j = 0; j < 12; j++) {
    if (nbr_del1[j] != -1)
      fprintf(fid, "1 %g\t%g\n", pos[nbr_del1[j]].x, pos[nbr_del1[j]].y);
  }
  fprintf(fid, "\n\n");

  fprintf(fid, "2 %g\t%g\n", pos[del2].x, pos[del2].y);
  for (int j = 0; j < 12; j++) {
    if (nbr_del2[j] != -1)
      fprintf(fid, "2 %g\t%g\n", pos[nbr_del2[j]].x, pos[nbr_del2[j]].y);
  }
  fprintf(fid, "\n\n");

  fprintf(fid, "3 %g\t%g\n", pos[add1].x, pos[add1].y);
  for (int j = 0; j < 12; j++) {
    if (nbr_add1[j] != -1)
      fprintf(fid, "3 %g\t%g\n", pos[nbr_add1[j]].x, pos[nbr_add1[j]].y);
  }
  fprintf(fid, "\n\n");

  fprintf(fid, "4 %g\t%g\n", pos[add2].x, pos[add2].y);
  for (int j = 0; j < 12; j++) {
    if (nbr_add2[j] != -1)
      fprintf(fid, "4 %g\t%g\n", pos[nbr_add2[j]].x, pos[nbr_add2[j]].y);
  }
  fprintf(fid, "\n\n");

  fflush(fid);
  fclose(fid);
  return 0;
}

double determinant(Vec3d X1, Vec3d X2, Vec3d X3, double len) {
  Vec3d A = diff_pbc(X1, X2, len);
  Vec3d B = diff_pbc(X1, X3, len);

  double det = (A.x * B.y - A.y * B.x);
  return det;
}



/* void wDiag(FILE *fid, MBRANE_para mbrane, AFM_para afm, SPRING_para spring, MESH mesh, */
/*             int i, int num_moves, double *Et,Vec3d *afm_force, */
/*             Vec3d *spring_force, double vol_sph, Vec3d *Pos){ */
/*     fprintf(fid, " %d %d %g %g %g %g", i, num_moves, mbrane.tot_energy[0], Et[0], Et[1], Et[2]); */
/*     if(afm.icompute!=0){fprintf(fid, " %g", Et[3]);} */
/*     if (spring.icompute!=0){fprintf(fid, " %g", Et[5]);} */
/*     if(fabs(mbrane.coef_vol_expansion)>1e-16){fprintf(fid, " %g", Et[4]);} */
/*     if (fabs(mbrane.pressure)>1e-16){fprintf(fid, " %g", Et[6]);} */
/*     if (afm.icompute!=0){fprintf(fid, " %g %g %g", afm_force->x,afm_force->y,afm_force->z );} */
/*     if (spring.icompute!=0){fprintf(fid, " %g %g", spring_force[0].z,spring_force[1].z);} */
/*     fprintf(fid, " %g %g %g %g\n",vol_sph,Pos[mesh.nPole].z,Pos[mesh.sPole].z, */
/*       height_rms(Pos,mbrane)); */
/*     fflush(fid); */
/* } */
/*-----------------------------------------------*/
