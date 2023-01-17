#include "global.h"
#include "subroutine.h"
std::mt19937 rng2;

extern "C" void Membrane_listread(int *, double *, double *, double *, 
        double *, double *, double *, bool *, char *);

extern "C" void Ljbot_listread(double *, double *, double *, 
        double *,  char *);

extern "C" void  MC_listread(char *, double *, double *, bool *,
           int *, int *, char *);

extern "C"  void  Activity_listread(char *, double *, double *, char *);

extern "C"  void  afm_listread(double *, double *, double *, double *,
             char *);

extern "C"  void  spring_listread(int *, double *, double *, char *);



void init_system_random_pos(Vec2d *Pos,  double len, 
        int N, char *metric, int bdry_condt ){

    /// @brief Initializes the points on surface of sphere or flat plane
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param metric Topology of the surface "cart" for flat plane "sph" for
    /// sphere
    ///  @param len length of the domain;
    ///  @param N number of points;

    bool is_sph, is_cart;
    double dl;

    // this should be calculated or passed parameters
    int n_ghost;
    // remove it once debugged 

    is_sph = false;
    is_cart = false;
    if(strcmp(metric, "sph") == 0){
        is_sph = true;
    }
    if(strcmp(metric, "cart") == 0){
        is_cart = true;
    }

    if(!is_cart && !is_sph){
        fprintf(stderr, "Unknown metric, input should be: \n");
        fprintf(stderr, "a) cart for a 2d cartesian monte carlo \n");
        fprintf(stderr, "b) sph  for a monte carlo on the surface of a sphere\n");
        exit(0);
    }


     
    if(is_cart){
        switch (bdry_condt) {
            case 0:
            // This is a channel ; read bdry_condt as 0 
            // n_ghost has to be even for logic to work;
            n_ghost = 2*(int) sqrt(N);
            dl = (2*len/n_ghost);
            for(int i=0; i<n_ghost/2; i++){
                    Pos[i].x = i*dl;
                    Pos[i].y = 0.0; 
            }
            for(int i=n_ghost/2; i<n_ghost; i++){
                Pos[i].x = (i - n_ghost/2 + 0.5)*dl;
                /* Pos[i].x = (i - n_ghost/2)*(2*(len + 0.5)/n_ghost);  bp*/
                Pos[i].y = len; 
            }
            for(int i=n_ghost; i<N; i++){
                Pos[i].x = drand48()*len;
                Pos[i].y = drand48()*len;
            }
            break;
      case 1:
            // This is a frame ;
            // n_ghost has to be even for logic to work;
            n_ghost = 4*(int) sqrt(N);
            dl = (4*len/n_ghost);
            for(int i=0; i<n_ghost/4; i++){
                Pos[i].x = i*dl;
                Pos[i].y = 0.0; 
            }
            for(int i=n_ghost/4; i<n_ghost/2; i++){
                Pos[i].x = (i - n_ghost/4 + 0.5)*dl;
                Pos[i].y = len; 
            }
            for(int i=n_ghost/2; i<3*n_ghost/4; i++){
                Pos[i].x = 0.0;
                Pos[i].y = (i - n_ghost/2 + 0.5)*dl; 
            }

            for(int i=3*n_ghost/4; i<n_ghost; i++){
                Pos[i].x = len;
                Pos[i].y = (i - 3*n_ghost/4)*dl; 
                }
            for(int i=n_ghost; i<N; i++){
                Pos[i].x = drand48()*len;
                Pos[i].y = drand48()*len;
            }

            break;
       
            default:
                for(int i=0; i<N; i++){
                    Pos[i].x = drand48()*len;
                    Pos[i].y = drand48()*len;
                }
                Pos[0].x = drand48()*len;
                Pos[0].y = drand48()*len;
        }
    }
    if(is_sph){
        Pos[0].x = 0;
        Pos[0].y = 0;
        Pos[1].x = pi;
        Pos[1].y = 0;
        for(int i=2; i<N; i++){
            Pos[i].x = acos(2*drand48() - 1); 
            Pos[i].y = 2*pi*drand48();
        }
        Pos[2].x = acos(2*drand48() - 1); 
        Pos[2].y = 2*pi*drand48();

    }
    /* for(int i=0; i<N; i++){ */
    /*     printf("%lf %lf \n", Pos[i].x); */
    /*     Pos[i].x = drand48()*len; */
    /*     Pos[i].y = drand48()*len; */
    /* } */

}

void init_eval_lij_t0(Vec3d *Pos, MESH mesh, double *lij_t0,
         MBRANE_para *para, SPRING_para *spring){
    /// @brief evaluates distance between neighbouring points and stores in lij_t0
    ///  @param Pos array containing co-ordinates of all the particles
   /// @param lij_t0 initial distance between points of membrane
    ///  @param mesh mesh related parameters -- connections and neighbours
    /// information; 
    ///  @param para membrane related parameters 
    Vec3d dr;
    int i,j,k;
    int num_nbr, cm_idx, npairs;
    double sum_lij=0;
    double r0;
    npairs = 0;
    for(i = 0; i < para->N; i++){
        num_nbr = mesh.numnbr[i];
        cm_idx = mesh.nghst * i;
        for(k = cm_idx; k < cm_idx + num_nbr; k++) {
            j = mesh.node_nbr_list[k];
            dr = diff_pbc(Pos[i] , Pos[j], para->len);
            /* dr = Pos[j] - Pos[i]; */
            lij_t0[k] = sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
            sum_lij += sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
            npairs++;
            /* printf("%g %g %g %g %g \n", Pos[i].x, Pos[j].x, Pos[i].y, Pos[j].y, lij_t0[k]); */
        }
    }

    para->av_bond_len = sum_lij/npairs;
    r0=para->av_bond_len;
    spring->constant=para->coef_bend/(r0*r0);
    if(para->is_fluid){
        for(i = 0; i < mesh.nghst*para->N; i++) lij_t0[i] = para->av_bond_len;
    }
}

void init_read_parameters( MBRANE_para *mbrane, 
        AFM_para *afm, MCpara *mcpara, ActivePara *activity,
        SPRING_para *spring, string para_file){
    /// @brief read parameters from para_file 
    ///  @param mesh mesh related parameters -- connections and neighbours
    /// information; 
    ///  @param mbrane membrane related parameters
    ///  @param afm AFM related parameters
    ///  @param mcpara monte-carlo related parameters
    //
    //
    char temp_algo[char_len];
    char which_act[char_len], tmp_fname[char_len];

    sprintf(tmp_fname, "%s", para_file.c_str() );

    Membrane_listread(&mbrane->N, &mbrane->coef_bend,
            &mbrane->YY, &mbrane->coef_vol_expansion,
            &mbrane->sp_curv,  &mbrane->pressure, &mbrane->radius,
            &mbrane->is_fluid, tmp_fname);

    sprintf(tmp_fname, "%s", para_file.c_str() );
    Ljbot_listread(&mbrane->pos_bot_wall, &mbrane->sigma,
            &mbrane->epsilon, &mbrane->theta,
            tmp_fname);


    sprintf(tmp_fname, "%s", para_file.c_str() );
    spring_listread(&spring->icompute, &spring->nPole_eq_z, &spring->sPole_eq_z, tmp_fname);

    sprintf(tmp_fname, "%s", para_file.c_str() );
    afm_listread(&afm->tip_rad, &afm->tip_pos_z, &afm->sigma, &afm->epsilon,
             tmp_fname);


    sprintf(tmp_fname, "%s", para_file.c_str() );
    MC_listread(temp_algo, &mcpara->dfac, &mcpara->kBT, &mcpara->is_restart,
            &mcpara->tot_mc_iter, &mcpara->dump_skip, tmp_fname);

    mcpara->algo=temp_algo;

    sprintf(tmp_fname, "%s", para_file.c_str() );
    Activity_listread(which_act, &activity->minA, &activity->maxA, tmp_fname);

    activity->act = which_act;

   /* printf("filename = %s %s\n", para_file.c_str(),  tmp_fname); */

   mbrane->num_triangles = 2*mbrane->N - 4;
   mbrane->num_nbr = 3*mbrane->num_triangles;
   // mbrane->av_bond_len = sqrt(8*pi/(2*mbrane->N-4));
   // define the monte carlo parameters
   mcpara->one_mc_iter = 2*mbrane->N;
   mcpara->delta = sqrt(8*pi/(2*mbrane->N-4));
}
//
void write_param(string fname, MBRANE_para mbrane, MCpara mcpara, SPRING_para spring){
    double FvK = mbrane.YY*mbrane.radius*mbrane.radius/mbrane.coef_bend;
    ofstream paramfile;
    paramfile.open( fname );
    paramfile << "# =========== Model Parameters ==========" << endl
            << "# Foppl von Karman number: FvK = " << FvK << endl
            << "# Elasto-thermal number: ET = " << mcpara.kBT/mbrane.coef_bend*sqrt(FvK) << endl
            << "# average bond length: r0 = " << mbrane.av_bond_len << endl;
    if (spring.icompute!=0){
       paramfile << "# Spring constant: Ki = " << spring.constant << endl;    
    }
    paramfile.close();
}


void init_activity(ActivePara activity, int N){
    int i;
    std::uniform_real_distribution<> rand_real(activity.minA, activity.maxA);
    if(activity.act == "random"){
        for(i=0;i<N;i++) activity.activity[i] = rand_real(rng2);
    }
    if(activity.act == "constant"){
        for(i=0;i<N;i++) activity.activity[i] = activity.maxA;
    }  
}
