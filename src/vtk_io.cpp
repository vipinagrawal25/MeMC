#include<math.h>
#include <iostream>
#include <string.h>
#include <errno.h>
#include <fstream>
using namespace std;
/*******************************************/
void vtk_io(double *points, 
        int *triangles, 
        int Np, string filename, string dataname){
    char first_headers[] = "# vtk DataFile Version 2.0 \n grid, mydata\n ASCII \n DATASET POLYDATA \n";
    int num_triangles, i;
    ofstream fid(filename,ofstream::out);
    num_triangles = 2*Np - 4;
    fid << first_headers;
    fid << "POINTS\t" << Np << "\t" << "float" << endl;
    for(i = 0; i< 3*Np; i=i+3){
        fid << points[i] << "\t" << points[i+1] << "\t" << points[i+2] << endl;
    }
    fid << "POLYGONS\t"<< num_triangles << "\t" << 4*num_triangles << endl;
    for(i = 0; i< 3*num_triangles; i=i+3){
        fid << 3 << "\t" << triangles[i] << "\t" << triangles[i+1] << "\t" << triangles[i+2] << endl;
    }
    fid.close();
}

 void vtk_io_point_data(double *data, 
         int Np, string filename, 
         string dataname){
     ofstream fid(filename,ofstream::app);
     int  i;

     fid << "POINT_DATA\t" << Np << endl;
     fid << "SCALARS\t" << dataname <<"\t" << "float 1" << endl;
     fid << "LOOKUP_TABLE default \n";
     for(i = 0; i< Np; i=i+1){
             fid << data[i] << endl;
         }
     fid.close();
 }

 void vtk_io_cell_data(double *data, 
         int Np, string filename, 
         string dataname){

     ofstream fid(filename,ofstream::app);
     int  i;

     fid << "CELL_DATA\t" << Np << endl;
     fid << "SCALARS\t" << dataname <<"\t" << "float 1" << endl;
     fid << "LOOKUP_TABLE default \n";
 
     for(i = 0; i< Np; i=i+1){
         fid << data[i] << endl;
     }
     fid.close();
 }

