#include <string>
#include "Vector.h"
#include <iomanip>
#include <sstream>
#ifndef FILE_MISC_SEEN
#define FILE_MISC_SEEN
using namespace std;
/* -----------------------------------------------*/
template<typename T>  // This type of function definition can take any variable type. 
T absolute(T value);       // Quite interesting, isn't it ?
double SqEr(double Arr1[], double Arr2[],int ndim);
bool FileExists(const std::string &s);
void write_param(string fname);
void print(double*, int);
void print(int*, int);
void print(double *arr, int np, int pp);
void print(double*, int start, int skip, int end);
void add(double *added, double *arr1, double *arr2, int nn);
void substract(double *subsed, double *arr1, double *arr2, int nn);
double norm(double *y, int nn );
double norm(const double *y, int nn );
void downScale(double *yscaled, double *y, int factor, int nn );
void zeros(double *yzero, int ndim);
void max(int *amaxind, double *amaxval, Vec3d *pos, int ndim, char dirn = 'z');
void min(int *aminind, double *aminval, Vec3d *pos, int ndim,char dirn='z');
int get_nstart(int , int );
void print(double *arr, int nn, string fname);
/* void wHeader(FILE *fid, MBRANE_p mbrane, AFM_p afm, SPRING_p spring); */
/* void wDiag(FILE *fid, MBRANE_p mbrane, AFM_p afm, SPRING_p spring, MESH mesh, */
/*             int i, int num_moves, double *Et,Vec3d *afm_force, */
/*             Vec3d *spring_force, double vol_sph,Vec3d *Pos); */
/*-----------------------------------------------*/
template<typename T>
inline string ZeroPadNumber(T num){
    ostringstream ss;
    ss << setw( 5 ) << setfill( '0' ) << (int)num;
    return ss.str();
}
/********************************************************/
#endif
