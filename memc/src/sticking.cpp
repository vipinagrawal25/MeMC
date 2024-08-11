#include "sticking.hpp"
#include <fstream>

extern "C" void StickRead(double *, double *, double *, char *);

int get_nstart(int, int);

int STICK::initSTICK(int N, std::string fname){
  char tmp_fname[128];
  string parafile, outfile;

  parafile = fname+"/para_file.in";
  sprintf(tmp_fname, "%s", parafile.c_str() );
  StickRead(&epsilon, &sigma, &pos_bot_wall,  tmp_fname);

  ofstream out_;
  out_.open( fname+"/stickpara.out");
  out_<< "# =========== stretching parameters ==========" << endl
      << " N " << N << endl
      << " epsilon = " << epsilon << endl
      << " sigma " << sigma << endl
      << " pos_bot_wall " << pos_bot_wall << endl;
  out_.close();
  return 0;

}

double lj_attr(double sqdr, double eps){

    /// @param sqdr square of the distance between two points
    /// @param eps coefficient of the potential 
    /// @return Energy evaluated using Lennard-Jones potential.
    ///  @details see https://en.wikipedia.org/wiki/Lennard-Jones_potential

    double r6;
    r6 = sqdr*sqdr*sqdr;
    return 4*eps*(r6*(r6-1));
}


double STICK::stick_energy_ipart(Vec3d pos){
    /// @brief Sticking of the bottom surface using LJ potential 
    /// @param zz z-coordinate of a point in shell
    /// @param eps coefficient of the potential 
    /// @param sigma sigma in the expression of LJ potential 
    /// @param sur_pos position of the bottom wall 
    /// @param is_attractive true if the zz sees the  bottom wall 
    /// @return Energy evaluated using Lennard-Jones potential.
    ///  @details see https://en.wikipedia.org/wiki/Lennard-Jones_potential

    double inv_sqdz, ds;

        ds = pos_bot_wall - pos.z;
        inv_sqdz = (sigma*sigma)/(ds*ds);
        return  lj_attr(inv_sqdz, epsilon);
}



double STICK::stick_energy_total(Vec3d *pos, 
         int N){

    /// @brief Estimate the total Sticking  
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param is_attractive true for all the particles which sees bottom wall 
    ///  @param para  Membrane related parameters;
    /// @return Total Energy contribution from Sticking to bottom surface

    int idx;
    double lj_bote;

    lj_bote = 0e0;
    for(idx = 0; idx < N; idx++){
        /* idx = 2; */
        lj_bote += stick_energy_ipart(pos[idx]);
    }
    return lj_bote;
}
