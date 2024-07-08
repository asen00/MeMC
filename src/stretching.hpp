#ifndef STRETCHING_HPP
#define STRETCHING_HPP
#include <string>
#include <vector>
#include "mesh.hpp"
#include "vector.hpp"

class STE {
public : 
    double stretch_energy_total(Vec3d *pos, MESH_p mesh);
    double area_ipart(Vec3d *pos, double *area, int *node_nbr, int num_nbr, int idx);
    double init_eval_lij_t0(Vec3d *Pos, MESH_p mesh,  bool is_fluid);
    double stretch_energy_ipart(Vec3d *pos, int *node_nbr, int num_nbr, int idx, int );
    double area_ipart(Vec3d *pos, int *node_nbr, int num_nbr, int idx);
    double area_total(Vec3d *, MESH_p );
    double volume_total(Vec3d *, MESH_p );
    double volume_ipart(Vec3d *, int *, int , int );
    int initSTE(int , std::string );
    bool dovol(){return do_volume;}
    bool dopressure(){return is_pressurized;}
    double getpressure(){return pressure;}
    double PV_change(double dvol){return pressure*dvol;}
private:
    double YY;  //coefficient stretching
    bool do_volume;
    bool is_pressurized;
    double coef_vol_expansion;   //coefficient of volume expansion
    double pressure;
    bool do_area;
    double coef_area_expansion;   //coefficient of area expansion
    vector <double> area_t0;    //  area of each triangle at t=0
    vector <double> lij_t0;
};
#endif