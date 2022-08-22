#include "global.h"
#include "subroutine.h"
#include <random>
#include <unistd.h>
#include <cmath>
#include <sstream>
#include <iomanip>
std::mt19937 rng;

void init_rng(uint32_t seed_val){

	 ///  @brief Generates random number 
    rng.seed(seed_val);
}

bool Metropolis(double DE, double kbt ){
    /// @brief Metropolis algorithm
    /// @param DE change in energy
    /// @param kbt boltzmann constant times temperature
    /// @return True if DE< 0 or the random number generated is less than
    /// exp(-DE/kbt)
    /// @details see https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm 
    bool yes;
    double rand;
    std::uniform_real_distribution<> rand_real(0, 1);

    yes = (DE <= 0.e0);
    if (!yes){
        rand = rand_real(rng);
        yes = rand < exp(-DE/kbt);
    }
    return yes;
}

double rand_inc_theta(double th0, 
        double dfac){
    /// @brief increment the polar angle randomly  
    double dth;
    double tmp_th0;
    std::uniform_real_distribution<> rand_real(-1, 1);


    tmp_th0 = 10;
    while (tmp_th0 > pi || tmp_th0 < 0){
        dth = (pi/dfac)*(rand_real(rng)); 
        tmp_th0 = th0 + dth;
    }
    return dth;
}

double energy_mc_3d(Vec3d *pos, MESH mesh, 
        double *lij_t0, bool *is_attractive, int idx,
        MBRANE_para mbrane, 
        MCpara mcpara, AFM_para afm){

    /// @brief Estimate the contribution from all the energies when a particle is moved randomly 
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param mesh mesh related parameters -- connections and neighbours
    ///information
    /// @param lij_t0 initial distance between points of membrane
    /// @param is_attractive true if the zz sees the  bottom wall 
    ///  @param idx index of ith particle which is moved;
    ///  @param mbrane  Membrane related parameters;
    /// @param mcpara Monte-Carlo related parameters
    /// @param AFM afm related parameter 
    /// @return Change in Energy when idx particle is moved


    double E_b, E_s;
    double E_stick;
    double  E_afm;
    int cm_idx, num_nbr;

    E_b = 0; E_s = 0; E_stick = 0; E_afm = 0;
    num_nbr = mesh.cmlist[idx + 1] - mesh.cmlist[idx];
    cm_idx = mesh.cmlist[idx];
    E_b = bending_energy_ipart(pos, 
            (int *) (mesh.node_nbr_list + cm_idx),
             num_nbr, idx, mbrane);

    E_b += bending_energy_ipart_neighbour(pos, 
            mesh, idx, mbrane);

    E_s = stretch_energy_ipart(pos, 
            (int *) (mesh.node_nbr_list + cm_idx),
            (double *) (lij_t0 + cm_idx), num_nbr, 
            idx, mbrane);

    E_stick = lj_bottom_surface(pos[idx].z, is_attractive[idx], 
            mbrane.pos_bot_wall, mbrane.epsilon, mbrane.sigma);

    E_afm = lj_afm(pos[idx], afm);

    return E_b + E_s + E_stick + E_afm;
}

int monte_carlo_3d(Vec3d *pos, MESH mesh, 
                double *lij_t0, bool *is_attractive, 
                MBRANE_para mbrane, 
                MCpara mcpara, AFM_para afm){

    /// @brief Monte-Carlo routine for the membrane 
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param mesh mesh related parameters -- connections and neighbours
    ///information
    /// @param lij_t0 initial distance between points of membrane
    /// @param is_attractive true if the zz sees the  bottom wall 
    ///  @param mbrane  Membrane related parameters;
    /// @param mcpara Monte-Carlo related parameters
    /// @param AFM afm related parameter 
    /// @return number of accepted moves 

    int i, move;
    int num_nbr, cm_idx;
    double x_o, y_o, z_o, x_n, y_n, z_n;
    double de,  Eini, Efin;
    double dxinc, dyinc, dzinc;
    double vol_i, vol_f;
    double dvol, de_vol, ini_vol;
    double KAPPA;
    std::uniform_int_distribution<uint32_t> rand_int(0,mbrane.N-1);
    std::uniform_real_distribution<> rand_real(-1, 1);
    ini_vol = (4./3.)*pi*pow(mbrane.radius,3);
    KAPPA = mbrane.coef_vol_expansion;
    move = 0;
    for(i = 0; i< mcpara.one_mc_iter; i++){
        int idx = rand_int(rng);
        num_nbr = mesh.cmlist[idx + 1] - mesh.cmlist[idx];
        cm_idx = mesh.cmlist[idx];

        Eini = energy_mc_3d(pos, mesh, 
                lij_t0, is_attractive, 
                idx,  mbrane, mcpara, afm);

        vol_i = volume_ipart(pos, 
                (int *) (mesh.node_nbr_list + cm_idx),
                num_nbr, idx, mbrane);

        x_o =   pos[idx].x;
        y_o =   pos[idx].y;
        z_o =   pos[idx].z;

        dxinc = (mcpara.delta/mcpara.dfac)*(rand_real(rng));
        dyinc = (mcpara.delta/mcpara.dfac)*(rand_real(rng));
        dzinc = (mcpara.delta/mcpara.dfac)*(rand_real(rng));

        x_n = x_o + dxinc;
        y_n = y_o + dyinc;
        z_n = z_o + dzinc;

        pos[idx].x = x_n;
        pos[idx].y = y_n;
        pos[idx].z = z_n;

        Efin = energy_mc_3d(pos, mesh, 
                lij_t0, is_attractive, 
                idx,   mbrane, mcpara, afm);

        vol_f = volume_ipart(pos, 
                (int *) (mesh.node_nbr_list + cm_idx),
                num_nbr, idx, mbrane);

        dvol =  (vol_f - vol_i);
        de_vol = (2*dvol/(ini_vol*ini_vol))*(mbrane.volume[0]  - ini_vol)
            + (dvol/ini_vol)*(dvol/ini_vol);
        de_vol = KAPPA*de_vol;
        de = (Efin - Eini) + de_vol;
        if (Metropolis(de,mcpara.kBT)){
            move = move + 1;
            mbrane.tot_energy[0] +=  de;
            mbrane.volume[0] += dvol;
        }
        else{
            pos[idx].x = x_o;
            pos[idx].y = y_o;
            pos[idx].z = z_o;
        }
    }
    return move;
}

int monte_carlo_surf2d(Vec2d *Pos, 
        Nbh_list *neib, LJpara para, 
        MCpara mcpara, char *metric){

    /// @brief Monte-Carlo routine for the points on the surface of sphere or flat
    /// plane
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param neib Neighbour list information for the particle 
    ///  @param para parameters related to LJ potential 
    /// @param mcpara Monte-Carlo related parameters
    ///  @param metric Topology of the surface "cart" for flat plane "sph" for
    /// sphere
    /// @return number of accepted moves 


    int i, move;
    double x_o, y_o, x_n, y_n;
    double de,  Eini, Efin;
    double dxinc, dyinc;
    bool is_sph, is_cart;
    int bdry_condt, n_ghost;

    bdry_condt = 0;
    switch (bdry_condt){
        case 0:
            //This is the channel
            n_ghost = 2*(int)sqrt(para.N);
            break;
        default:
            n_ghost = 0;
    }
    std::uniform_int_distribution<uint32_t> rand_int(n_ghost,para.N-1);
    std::uniform_real_distribution<> rand_real(-1, 1);

    is_sph = false;
    is_cart = false;

    if(strcmp(metric, "sph") == 0){
        is_sph = true;
    }
    if(strcmp(metric, "cart") == 0){
        is_cart = true;
    }


    move = 0;

    for(i = 0; i< mcpara.one_mc_iter; i++){
        int idx = rand_int(rng);
        Eini =  pairlj_ipart_energy(Pos, neib[idx].list,
                neib[idx].cnt, idx, para, metric);
        x_o =   Pos[idx].x;
        y_o =   Pos[idx].y;

        if(is_cart){
            dxinc = (para.sigma/mcpara.dfac)*(2*rand_real(rng) - 1);
            dyinc = (para.sigma/mcpara.dfac)*(2*rand_real(rng) - 1);
            x_n = fmod((x_o + dxinc + 30*para.len), para.len);
            y_n = fmod((y_o + dyinc + 30*para.len), para.len);
            Pos[idx].x = x_n;
            Pos[idx].y = y_n;
        }
        if(is_sph){
            dxinc = rand_inc_theta(Pos[idx].x, mcpara.dfac);
            dyinc = (para.sigma/mcpara.dfac)*(2*rand_real(rng) - 1);
            x_n = x_o + para.sigma*dxinc;
            y_n = fmod((y_o + dyinc + 30*2*pi), 2*pi);
            Pos[idx].x = x_n;
            Pos[idx].y = y_n;
        }

        Efin =  pairlj_ipart_energy(Pos, neib[idx].list,
                neib[idx].cnt, idx, para, metric);
        de = (Efin - Eini);
        if(Metropolis(de, mcpara.kBT)){
            move = move + 1;
        }
        else{
            Pos[idx].x = x_o;
            Pos[idx].y = y_o;
        } 
    }
    return move;
}

