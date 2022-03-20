#include "global.h"
#include "subroutine.h"
#include <random>
#include <unistd.h>
#include <cmath>
#include <sstream>
#include <iomanip>
std::mt19937 rng;

void init_rng(uint32_t seed_val){
    rng.seed(seed_val);
}

bool Metropolis(double DE, double kbt ){
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
        double *lij_t0, bool *is_attractive, int idx, bool *is_be_pos,
        MBRANE_para mbrane, 
        MCpara mcpara, AFM_para afm){
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
    *is_be_pos = E_b > 0;

    return E_b + E_s + E_stick + E_afm;
}

int monte_carlo_3d(Vec3d *pos, MESH mesh, 
                double *lij_t0, bool *is_attractive, 
                MBRANE_para mbrane, 
                MCpara mcpara, AFM_para afm){
    int i, move;
    int num_nbr, cm_idx;
    double x_o, y_o, z_o, x_n, y_n, z_n;
    double de,  Eini, Efin;
    double dxinc, dyinc, dzinc;
    double vol_i, vol_f;
    double dvol, de_vol, ini_vol;
    double KAPPA;
    bool is_be_pos;
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
                idx, &is_be_pos, mbrane, mcpara, afm);

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
                idx, &is_be_pos,  mbrane, mcpara, afm);

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
    int i, move;
    double x_o, y_o, x_n, y_n;
    double de,  Eini, Efin;
    double dxinc, dyinc;
    bool is_sph, is_cart;

    std::uniform_int_distribution<uint32_t> rand_int(0,para.N-1);
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

