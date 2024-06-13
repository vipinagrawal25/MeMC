#include<math.h>
#include <iostream>
#include <string.h>
#include <errno.h>
#include <fstream>
using namespace std;
/*******************************************/
void visit_vtk_io(double *points, 
        int *triangles, 
        int Np, string filename, string dataname){
    char first_headers[] = "# vtk DataFile Version 2.0 \n grid, mydata\n ASCII \n DATASET POLYDATA \n";
    int num_triangles, i;
    ofstream fid(filename,ofstream::out);
    // FILE *fid;
    // //
    num_triangles = 2*Np - 4;
    // fid = fopen(filename, "w");
    // if (fid==NULL){
    //     cout <<"FILE ERROR: can not open the " << filename << endl;
    //     cout <<"Error = " << errno << endl;
    // }
    fid << first_headers;
    fid << "POINTS\t" << Np << "\t" << "float" << endl;
    // fprintf(fid, "%s", first_headers);
    // fprintf(fid, "%s %d %s", "POINTS", Np, "float \n");
    for(i = 0; i< 3*Np; i=i+3){
        fid << points[i] << "\t" << points[i+1] << "\t" << points[i+2] << endl;
        // fprintf(fid, "%g %g %g \n", points[i], points[i+1], points[i+2]);
    }
    fid << "POLYGONS\t"<< num_triangles << "\t" << 4*num_triangles << endl;
    // fprintf(fid, "%s %d  %d %s", "POLYGONS", num_triangles, 4*num_triangles, " \n");
    for(i = 0; i< 3*num_triangles; i=i+3){
        fid << 3 << "\t" << triangles[i] << "\t" << triangles[i+1] << "\t" << triangles[i+2] << endl;
        // fprintf(fid, "%d %d %d %d\n", 3, triangles[i], triangles[i+1], triangles[i+2]);
    }
    // fclose(fid);
    fid.close();
}
void visit_vtk_io(double *points, 
        int Np, string filename){
    '''
        Dumps data into .vtk file
    '''
    char first_headers[] = "# vtk DataFile Version 2.0 \n grid, mydata\n ASCII \n DATASET POLYDATA \n";
    int i;
    ifstream fid(filename,ofstream::out);
    fid >> first_headers;
    fid << "POINTS\t" << Np << "\t" << "float" << endl;
    for(i = 0; i< 3*Np; i=i+3){
        fid << points[i] << "\t" << points[i+1] << "\t" << points[i+2] << endl;
    }
    fid.close();
}
void visit_vtk_read(double *points, int Np, string filename){
    ''' 
        Function to read data from .vtk files.
    '''
    ifstream myfile;
    myfile.open(fname);
    double ch;
    while(myfile>>ch)
}

void visit_vtk_io_point_data(bool *data, 
        int Np, string filename, 
        string dataname){
    ofstream fid(filename,ofstream::app);
    int num_triangles, i;
    num_triangles = 2*Np - 4;
    fid << "POINT_DATA\t" << Np << endl;
    fid << "SCALARS\t" << dataname <<"\t" << "float 1" << endl;
    fid << "LOOKUP_TABLE default \n";
    // fprintf(fid, "%s %d \n", "POINT_DATA", Np);
    // fprintf(fid, "%s %s %s \n", "SCALARS", dataname, "float 1");
    // fprintf(fid, "%s \n", "LOOKUP_TABLE default ");
    for(i = 0; i< Np; i=i+1){
        if(data[i]){
            fid << 1.0 << endl;
            // fprintf(fid, "%g \n", 1.0);
        } else {
            fid << 0.0 << endl;
            // fprintf(fid, "%g \n", 0.0);
        }
    }
    fid.close();
}

// void visit_vtk_io_cell_data(double *data, 
//         int Np, string filename, 
//         string dataname){

//     int num_triangles, i;
//     FILE *fid;

//     num_triangles = 2*Np - 4;
//     fid = fopen(filename, "a");
//     /* if(fid != NULL)printf("I am null"); */
//     fprintf(fid, "%s %d  \n", "CELL_DATA", Np);
//     fprintf(fid, "%s %s %s \n", "SCALARS", dataname, "float 1");
//     fprintf(fid, "%s \n", "LOOKUP_TABLE default ");
//     for(i = 0; i< Np; i=i+1){
//         fprintf(fid, "%g \n", data[i]);
//     }
//     fclose(fid);
// }

// void visit_vtk_io_afm_tip(double *data, 
//         int Np, string filename){

//     int i, j, here, est, nth, n_est;
//     double x, y;
//     FILE *fid;
//     int iext, idx;
//     double dx;

//     iext = (int) sqrt(Np);
//     char first_headers[] = "# vtk DataFile Version 2.0 \n grid, afmtip\n ASCII \n DATASET POLYDATA \n";

//     fid = fopen(filename, "w");
//     fprintf(fid, "%s", first_headers);
//     fprintf(fid, "%s %d %s", "POINTS", Np, "float \n");

//     for(i=0; i < 3*Np; i=i+3){
//         fprintf(fid, "%g %g %g\n", data[i], data[i+1], data[i+2]);
//     }

//     int npoly = (iext-1)*(iext-1);
//     fprintf(fid, "%s %d  %d %s", "POLYGONS", 4*npoly, 5*4*npoly, " \n");
//     for(j=0; j < iext-1; j++){
//         for(i=0; i < iext-1; i++){
//             here = j*iext + i;
//             est = (i+1)%iext + ((iext + j)%iext) *iext;
//             nth = ((iext + i)%iext) + (j + 1 + iext)%iext *iext;
//             n_est = ((i + 1)%iext) + ((j + 1 + iext)%iext)*iext;
//             fprintf(fid, "%d %d %d %d %d\n", 4, here, est, n_est, nth);
//         }
//     }
//     fclose(fid);
// }