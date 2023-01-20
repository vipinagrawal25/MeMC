#include "global.h"
#include "subroutine.h"
#include <cstdio>
#include <random>
#include <unistd.h>
#include <cmath>
#include <sstream>
#include <iomanip>
std::mt19937 rng;

int print_sanity(int *nbr_del1, int *nbr_del2, 
        int *nbr_add1, int *nbr_add2, int del1, int del2,
        int add1, int add2, char *fname){
    FILE *fid;
    int k=10;

    fid = fopen(fname, "a");
    /* fprintf(fid, "start\n"); */
    fprintf(fid, "%06d\t%06d\t%06d\t%06d\t%06d\n", k, del1, del2, add1, add2); 
    for(int j=0; j<12; j++){
        fprintf(fid, "%06d\t%06d\t%06d\t%06d\t%06d\n",j, nbr_del1[j],
                nbr_del2[j], nbr_add1[j], nbr_add2[j]);
    }
    fprintf(fid, "next line\n"); 
    fflush(fid);
    fclose(fid);
    return 0;
}


int print_sanity(Vec3d *pos, int *nbr_del1, int *nbr_del2, 
        int *nbr_add1, int *nbr_add2, int del1, int del2,
        int add1, int add2, char *ftag, int idx){
    FILE *fid;
    char fname[128];
    int k=10;

    sprintf(fname, "check_pos/%s_pos_%04d_.dat",ftag,idx);
    fid = fopen(fname, "w");
    /* fprintf(fid, "start\n"); */
    fprintf(fid, "1 %g\t%g\n", pos[del1].x, pos[del1].y); 
    for(int j=0; j<12; j++){
        if(nbr_del1[j] != -1)
        fprintf(fid, "1 %g\t%g\n",pos[nbr_del1[j]].x, pos[nbr_del1[j]].y);
    }
    fprintf(fid, "\n\n"); 

    fprintf(fid, "2 %g\t%g\n", pos[del2].x, pos[del2].y); 
    for(int j=0; j<12; j++){
        if(nbr_del2[j] != -1)
        fprintf(fid, "2 %g\t%g\n",pos[nbr_del2[j]].x, pos[nbr_del2[j]].y);
    }
    fprintf(fid, "\n\n"); 


    fprintf(fid, "3 %g\t%g\n", pos[add1].x, pos[add1].y); 
    for(int j=0; j<12; j++){
        if(nbr_add1[j] != -1)
        fprintf(fid, "3 %g\t%g\n",pos[nbr_add1[j]].x, pos[nbr_add1[j]].y);
    }
    fprintf(fid, "\n\n"); 


    fprintf(fid, "4 %g\t%g\n", pos[add2].x, pos[add2].y); 
    for(int j=0; j<12; j++){
        if(nbr_add2[j] != -1)
        fprintf(fid, "4 %g\t%g\n",pos[nbr_add2[j]].x, pos[nbr_add2[j]].y);
    }
    fprintf(fid, "\n\n"); 

    fflush(fid);
    fclose(fid);
    return 0;
}


double determinant(Vec3d X1, Vec3d X2, Vec3d X3, double len){
    Vec3d A = diff_pbc(X1, X2, len);
    Vec3d B = diff_pbc(X1, X3, len);

    double det = (A.x*B.y - A.y*B.x);
return det;
}


void init_rng(uint32_t seed_val){

	 ///  @brief Generates random number 
    rng.seed(seed_val);
}

int del_nbr(int *nbrs, int numnbr, int idx){
    // delet int idx between i1 and i2 in the nbrs list
    int new_numnbr, delete_here;
    bool logic;

    logic = false;
    delete_here = 0;

    // for(int i=0; i<numnbr+3; i++)printf("%d \n", nbrs[i]);
        // printf("\n\n");


    while(!logic){
            logic = (nbrs[delete_here] == idx);
            ++delete_here;
        }

    memcpy(nbrs+delete_here-1, &nbrs[delete_here], sizeof(int) * (numnbr-delete_here+1));

    // for(int i=0; i<numnbr+3; i++)printf("%d \n", nbrs[i]);
        // printf("\n\n");

    return numnbr-1;
}


int add_nbr(int *nbrs, int numnbr, int idx, int i1, int i2){
    // add int idx between i1 and i2 in the nbrs list
    int insert_here;
    bool logic;

    logic = false;
    insert_here = 0;

    // for(int i=0; i<numnbr+3; i++)printf("%d \n", nbrs[i]);
    //     printf("\n\n");

    while(!logic){
            logic = (nbrs[insert_here] == i1) || (nbrs[insert_here] == i2);
            ++insert_here;
        }

        logic = (nbrs[insert_here] == i1) || (nbrs[insert_here] == i2);
        if (logic){
            memcpy(nbrs+insert_here, &nbrs[insert_here-1], sizeof(int) * (numnbr-insert_here+1));
            nbrs[insert_here] = idx;
        }else{
            insert_here = 0;
            memcpy(nbrs+insert_here, &nbrs[insert_here-1], sizeof(int) * (numnbr-insert_here+1));
            nbrs[insert_here] = idx;
        }

    // for(int i=0; i<numnbr+3; i++)printf("%d \n", nbrs[i]);
    //     printf("\n\n");

    return numnbr+1;
}

bool Metropolis(double DE, double activity, MCpara mcpara){
    /// @brief Metropolis algorithm
    /// @param DE change in energy
    /// @param kbt boltzmann constant times temperature
    /// @return True if DE< 0 or the random number generated is less than
    /// exp(-DE/kbt)
    /// @details see https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm 
    bool yes;
    double rand;
    DE += activity;
    std::uniform_real_distribution<> rand_real(0, 1);
    yes = (DE <= 0.e0);
    if (!yes){
        rand = rand_real(rng);
        yes = rand < exp(-DE/mcpara.kBT);
    }
    return yes;
}
bool Glauber(double DE, double activity, MCpara mcpara){
    /// @brief Glauber algorithm
    /// @param DE change in energy
    /// @param kbt boltzmann constant times temperature
    bool yes;
    double rand;
    DE += activity;
    std::uniform_real_distribution<> rand_real(0, 1);
    rand = rand_real(rng);
    yes=rand < 1/(1+exp(DE/mcpara.kBT));
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
        MBRANE_para mbrane, MCpara mcpara, AFM_para afm,
         SPRING_para spring){
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


    double E_b, E_s, E_stick, E_afm, E_spr;
    int cm_idx, num_nbr;

    E_b = 0; E_s = 0; E_stick = 0; E_afm = 0;
    
    cm_idx = mesh.nghst*idx;
    num_nbr = mesh.numnbr[idx];
   
    E_b = bending_energy_ipart(pos, 
            (int *) (mesh.node_nbr_list + cm_idx),
             num_nbr, idx, mbrane);

    E_b += bending_energy_ipart_neighbour(pos, 
            mesh, idx, mbrane);

    E_s = stretch_energy_ipart(pos, 
            (int *) (mesh.node_nbr_list + cm_idx),
            (lij_t0 + cm_idx), num_nbr,
            idx, mbrane);

    /* E_stick = lj_bottom_surface(pos[idx].z, is_attractive[idx], */ 
            /* mbrane.pos_bot_wall, mbrane.epsilon, mbrane.sigma); */

    /* E_afm = lj_afm(pos[idx], afm); */

    /* E_spr = spring_energy(pos[idx], idx, mesh, spring); */
    return E_b + E_s; // + E_stick + E_afm + E_spr;
}

int monte_carlo_3d(Vec3d *pos, MESH mesh, 
                double *lij_t0, bool *is_attractive, 
                MBRANE_para mbrane, MCpara mcpara, AFM_para afm, 
                ActivePara activity, SPRING_para spring){

    /// @brief Monte-Carlo routine for the membrane 
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param mesh mesh related parameters -- connections and neighbours
    /// information
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
    double dvol, de_vol, ini_vol, de_pressure;
    double KAPPA;
    bool yes;
    
    int nframe;
    nframe = 4*(int)sqrt(mbrane.N);

    std::uniform_int_distribution<uint32_t> rand_int(nframe,mbrane.N-1);
    std::uniform_real_distribution<> rand_real(-1, 1);
    ini_vol = (4./3.)*pi*pow(mbrane.radius,3);
    KAPPA = mbrane.coef_vol_expansion;
    move = 0;

    for(i = 0; i< mcpara.one_mc_iter; i++){
        int idx = rand_int(rng);
        /*
        num_nbr = mesh.cmlist[idx + 1] - mesh.cmlist[idx];
        cm_idx = mesh.cmlist[idx];
        */

        Eini = energy_mc_3d(pos, mesh, 
                lij_t0, is_attractive, 
                idx,  mbrane, mcpara, afm, spring);

        /* vol_i = volume_ipart(pos, */ 
        /*         (int *) (mesh.node_nbr_list + cm_idx), */
        /*         num_nbr, idx, mbrane); */

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
                idx,   mbrane, mcpara, afm, spring);

        /* vol_f = volume_ipart(pos, */ 
        /*         (int *) (mesh.node_nbr_list + cm_idx), */
        /*         num_nbr, idx, mbrane); */

        /* dvol=0.5*(vol_f - vol_i); */
        /* de_vol = vol_energy_change(mbrane,dvol); */
        /* de_pressure = PV_change(mbrane,dvol); */
        // de_vol = (2*dvol/(ini_vol*ini_vol))*(mbrane.volume[0]  - ini_vol)
        //     + (dvol/ini_vol)*(dvol/ini_vol);
        // de_vol = KAPPA*de_vol;

        de = (Efin - Eini); //+ de_vol + de_pressure;
        if (mcpara.algo == "mpolis"){
            yes=Metropolis(de, activity.activity[idx], mcpara);
        }else if(mcpara.algo == "glauber"){
            yes=Glauber(de, activity.activity[idx], mcpara);
        }

        if (yes){
            move = move + 1;
            mbrane.tot_energy[0] +=  de;
            /* mbrane.volume[0] += dvol; */
        }else{
            pos[idx].x = x_o;
            pos[idx].y = y_o;
            pos[idx].z = z_o;
        }
    }
    return move;
}

int monte_carlo_surf2d(Vec2d *Pos, Nbh_list *neib, LJpara para, 
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
    int n_ghost;

    switch (para.bdry_condt){
        case 0:
            n_ghost = 2*(int)sqrt(para.N);
            break;
        case 1:
            n_ghost = 4*(int)sqrt(para.N);
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
            dxinc = (para.sigma/mcpara.dfac)*rand_real(rng);
            x_n = fmod((x_o + dxinc + 30*para.len), para.len);
            Pos[idx].x = x_n;

            dyinc = (para.sigma/mcpara.dfac)*rand_real(rng);
            y_n = fmod((y_o + dyinc + 30*para.len), para.len);
            Pos[idx].y = y_n;
        }
        if(is_sph){
            dxinc = rand_inc_theta(Pos[idx].x, mcpara.dfac);
            dyinc = (para.sigma/mcpara.dfac)*(rand_real(rng));
            x_n = x_o + para.sigma*dxinc;
            y_n = fmod((y_o + dyinc + 30*2*pi), 2*pi);
            Pos[idx].x = x_n;
            Pos[idx].y = y_n;
        }

        Efin =  pairlj_ipart_energy(Pos, neib[idx].list,
                neib[idx].cnt, idx, para, metric);
        de = (Efin - Eini);
        if(Metropolis(de, 0.0,  mcpara)){
            move = move + 1;
        }
        else{
            Pos[idx].x = x_o;
            Pos[idx].y = y_o;
        } 
    }
    return move;
}

int monte_carlo_fluid(Vec3d *pos, MESH mesh,
                MBRANE_para mbrane, MCpara mcpara, AFM_para afm,
                ActivePara activity, SPRING_para spring){

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

    int i, j, move;
    int nnbr_del1;
    int cm_idx_del1, cm_idx_del2;
    int cm_idx_add1, cm_idx_add2;
    int idx_del1, idx_del2;
    int idx_add1, idx_add2;
    int nframe;

    int nbr_add_1[12], nbr_add_2[12];
    int nbr_del_1[12], nbr_del_2[12];

    int N_nbr_del2, N_nbr_del1;
    int N_nbr_add2, N_nbr_add1;
    double det1, det2;
    Vec3d bef_ij, aft_ij;

    double KAPPA;
    bool yes, logic;

    nframe = 4*(int)sqrt(mbrane.N);

    std::uniform_int_distribution<uint32_t> rand_int(nframe,mbrane.N-1);
    std::uniform_int_distribution<uint32_t> rand_nbr(0,mesh.nghst-1);
    std::uniform_real_distribution<> rand_real(-1, 1);

    move = 0;

    int idxn, up, down;

    for(i = 0; i< mcpara.one_mc_iter; i++){
        // identify the pair to be divorced
        // stored as idx_del1 and idx_del2
        logic = false;
        while (!logic){
            idx_del1 = rand_int(rng);
            cm_idx_del1 = mesh.nghst*idx_del1;
            nnbr_del1 = mesh.numnbr[idx_del1];
            idxn = rand_nbr(rng);
           if(mesh.node_nbr_list[cm_idx_del1 + idxn] != -1){
               idx_del2 = mesh.node_nbr_list[cm_idx_del1 + idxn];
               cm_idx_del2 = mesh.nghst*idx_del2;
                up = (idxn+1+nnbr_del1)%nnbr_del1;
                down = (idxn-1+nnbr_del1)%nnbr_del1;
                idx_add1 = mesh.node_nbr_list[idx_del1*mesh.nghst + up];
                idx_add2 = mesh.node_nbr_list[idx_del1*mesh.nghst + down];
                logic = idx_del2 > nframe && idx_add1 > nframe && idx_add2 > nframe; 
            }
            else{
                logic = false;
            }
        }


        det1 = determinant(pos[idx_add1], pos[idx_add2], pos[idx_del1], mbrane.len);
        det2 = determinant(pos[idx_add1], pos[idx_add2], pos[idx_del2], mbrane.len);
        cm_idx_add1 = mesh.nghst*idx_add1;
        cm_idx_add2 = mesh.nghst*idx_add2;


        /* print_sanity(mesh.node_nbr_list+cm_idx_del1, mesh.node_nbr_list+cm_idx_del2, */
        /*         mesh.node_nbr_list+cm_idx_add1, mesh.node_nbr_list+cm_idx_add2, */
        /*         idx_del1, idx_del2, idx_add1, idx_add2, "check.dat"); */

        if(det1*det2 < 0.0){
                    // copy the neighbours
            memcpy(nbr_del_1, &mesh.node_nbr_list[cm_idx_del1], sizeof(int) * mesh.nghst);
            memcpy(nbr_del_2, &mesh.node_nbr_list[cm_idx_del2], sizeof(int) * mesh.nghst);
            memcpy(nbr_add_1, &mesh.node_nbr_list[cm_idx_add1], sizeof(int) * mesh.nghst);
            memcpy(nbr_add_2, &mesh.node_nbr_list[cm_idx_add2], sizeof(int) * mesh.nghst);

            int is_true = 0;
            
            // form the bond
            N_nbr_add1 =  add_nbr(nbr_add_1, mesh.numnbr[idx_add1], idx_add2, idx_del1, idx_del2);
            N_nbr_add2 = add_nbr(nbr_add_2, mesh.numnbr[idx_add2], idx_add1, idx_del1, idx_del2);

            // get divorced
            N_nbr_del1 = del_nbr(nbr_del_1, mesh.numnbr[idx_del1], idx_del2);
            N_nbr_del2 = del_nbr(nbr_del_2, mesh.numnbr[idx_del2], idx_del1);

           
            bef_ij = pos[idx_del2] - pos[idx_del1];
            aft_ij = pos[idx_add2] - pos[idx_add1];
            /* bef_ij = diff_pbc(pos[idx_del1] , pos[idx_del2], mbrane.len); */
            /* aft_ij = diff_pbc(pos[idx_add1] , pos[idx_add2], mbrane.len); */

            double de = (norm(aft_ij) - mbrane.av_bond_len)*(norm(aft_ij) - mbrane.av_bond_len) - 
                (norm(bef_ij) - mbrane.av_bond_len)*(norm(bef_ij) - mbrane.av_bond_len); 

            // de = (Efin - Eini) + de_vol + de_pressure;
            de = mbrane.YY*de;
            if (mcpara.algo == "mpolis"){
                yes=Metropolis(de, 0.0e0, mcpara);
            }else if(mcpara.algo == "glauber"){
                yes=Glauber(de, 0.0e0, mcpara);
            }
            /* yes = true; */
            if (yes){
                move = move + 1;
                mbrane.tot_energy[0] +=  de;

                /* print_sanity(pos, mesh.node_nbr_list+cm_idx_del1, mesh.node_nbr_list+cm_idx_del2, */
                /*         mesh.node_nbr_list+cm_idx_add1, mesh.node_nbr_list+cm_idx_add2, */
                /*         idx_del1, idx_del2, idx_add1, idx_add2, (char*)"bef", i); */


                memcpy(mesh.node_nbr_list+cm_idx_del1, &nbr_del_1, sizeof(int) * mesh.nghst);
                memcpy(mesh.node_nbr_list+cm_idx_del2, &nbr_del_2, sizeof(int) * mesh.nghst);
                mesh.numnbr[idx_del1] = N_nbr_del1;
                mesh.numnbr[idx_del2] = N_nbr_del2;

                memcpy(mesh.node_nbr_list+cm_idx_add1, &nbr_add_1, sizeof(int) * mesh.nghst);
                memcpy(mesh.node_nbr_list+cm_idx_add2, &nbr_add_2, sizeof(int) * mesh.nghst);
                mesh.numnbr[idx_add1] = N_nbr_add1;
                mesh.numnbr[idx_add2] = N_nbr_add2;

                /* print_sanity(pos, mesh.node_nbr_list+cm_idx_del1, mesh.node_nbr_list+cm_idx_del2, */
                /*         mesh.node_nbr_list+cm_idx_add1, mesh.node_nbr_list+cm_idx_add2, */
                /*         idx_del1, idx_del2, idx_add1, idx_add2, (char*)"aft", i); */


                /* exit(0); */

            }
        }
    }
    return move;
}
