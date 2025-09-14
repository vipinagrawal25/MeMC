#ifndef STRETCHING_HPP
#define STRETCHING_HPP
#include <string>
#include <vector>
#include "mesh.hpp"
#include "vector.hpp"

class STE {
public :
    STE(const MESH_p&, std::string);
    double stretch_energy_total(Vec3d *pos, MESH_p mesh);
    void init_eval_lij_t0(MESH_p &mesh,  bool is_fluid);
    double stretch_energy_ipart(Vec3d *pos, int *node_nbr, int num_nbr, int idx, int, double lenth, int edge, bool pbc);
    double stretch_energy_ipart(double *, int, int, int);
    double area_total(MESH_p);
    double volume_total(Vec3d *, MESH_p );
    double volume_ipart(Vec3d *, int *, int , int, double lenth,
            int edge, bool pbc );
    double vol_energy_change(double volume, double dvol);
    bool doarea(){return do_area;}
    bool dovol(){return do_volume;}
    bool dopressure(){return is_pressurized;}
    double getpressure(){return pressure;}
    double PV_change(double dvol){return pressure*dvol;}
    double getkappa(){return Kappa;}
    void init_coefstretch(MESH_p);
    double area_energy_ipart(Vec3d *pos, int *node_nbr, int num_nbr, int idx,
        double lenth, int edge, bool pbc);
    double area_energy_total(MESH_p mesh);
    double getyy1(){return YY1;}
    double getyy2(){return YY2;}
private:
    double YY1, YY2;                // coefficient stretching
    bool do_volume;
    bool is_pressurized;
    double Kappa;                   //coefficient of volume expansion
    double pressure;
    bool do_area;
    double coef_area_expansion;     //coefficient of area expansion
    double *area_t0;                // area of each triangle at t=0
    vector <double> lij_t0;
    vector <double> HH;
    double ini_vol;
    double area_ipart(Vec3d *pos, int *node_nbr, int num_nbr, 
                    int idx, double lenth, int edge, bool pbc);
    void area_ipart(double* area, Vec3d *pos, int *node_nbr, int num_nbr, 
                    int idx, double lenth, int edge, bool pbc);
    void init_area_t0(MESH_p mesh);
    string initial_l0;
};
#endif