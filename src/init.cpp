#include "global.h"
#include "subroutine.h"
std::mt19937 rng2;

void init_system_random_pos(Vec2d *Pos,  double len, int N, char *metric){

    /// @brief Initializes the points on surface of sphere or flat plane
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param metric Topology of the surface "cart" for flat plane "sph" for
    /// sphere
    ///  @param len length of the domain;
    ///  @param N number of points;

    bool is_sph, is_cart;

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
        for(int i=0; i<N; i++){
            Pos[i].x = drand48()*len;
            Pos[i].y = drand48()*len;
        }
        Pos[0].x = drand48()*len;
        Pos[0].y = drand48()*len;
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
}

void init_eval_lij_t0(Vec3d *Pos, MESH mesh, 
        double *lij_t0, MBRANE_para *para, SPRING_para *spring){
    /// @brief evaluates distance between neighbouring points and stores in lij_t0
    ///  @param Pos array containing co-ordinates of all the particles
   /// @param lij_t0 initial distance between points of membrane
    ///  @param mesh mesh related parameters -- connections and neighbours
    /// information; 
    ///  @param para membrane related parameters 
    Vec3d dr;
    int i,j,k;
    int num_nbr, cm_idx;
    double sum_lij=0;
    double r0;
    for(i = 0; i < para->N; i++){
        num_nbr = mesh.cmlist[i + 1] - mesh.cmlist[i];
        cm_idx = mesh.cmlist[i];
        for(k = cm_idx; k < cm_idx + num_nbr; k++) {
            j = mesh.node_nbr_list[k];
            dr.x = Pos[i].x - Pos[j].x;
            dr.y = Pos[i].y - Pos[j].y;
            dr.z = Pos[i].z - Pos[j].z;
            lij_t0[k] = sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
            sum_lij+=lij_t0[k];
        }
    }
    para->av_bond_len = sum_lij/mesh.cmlist[para->N];
    r0=para->av_bond_len;
    spring->constant=para->coef_bend/(r0*r0);
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
    char buff[255];
    char temp_algo[char_len];
    int t_n, t_n2, t_n3, err;
    double td1, td2, td3, td4, td5, td6;
    char which_act[char_len];
    FILE *f2;
    f2 = fopen(para_file.c_str(), "r");
    if(f2){
        if(fgets(buff,255,(FILE*)f2) != NULL); 
        if(fgets(buff,255,(FILE*)f2) != NULL); 
        if(fgets(buff,255,(FILE*)f2) != NULL){
            sscanf(buff,"%d %lf %lf %lf %lf %lf %lf", 
                &t_n, &td1, &td2, &td3, &td4, &td5, &td6);
        }
        /* fprintf(stderr, "%s\n", buff); */
        mbrane->N = t_n;
        mbrane->coef_bend = td1;
        mbrane->YY = td2;       // Young's modulus        
        mbrane->coef_vol_expansion = td3;
        mbrane->sp_curv = td4;
        mbrane->pressure = td5;
        mbrane->radius = td6;
        if( fgets(buff,255,(FILE*)f2) != NULL);
        if( fgets(buff,255,(FILE*)f2) != NULL);
        if( fgets(buff,255,(FILE*)f2) != NULL);{
            sscanf(buff,"%lf %lf %lf %lf", &td1,&td2,&td3,&td4);
        }
        /* fprintf(stderr, "%s\n", buff); */
        mbrane->pos_bot_wall = td1;
        mbrane->sigma = td2;
        mbrane->epsilon = td3;
        mbrane->theta = td4;
        if (fabs(mbrane->epsilon)<1e-16 || fabs(mbrane->theta) < 1e-16){
            mbrane->istick = 0;        // second line done   
        }
        if( fgets(buff,255,(FILE*)f2) != NULL);
        if( fgets(buff,255,(FILE*)f2) != NULL); 
        if( fgets(buff,255,(FILE*)f2) != NULL){ 
            sscanf(buff,"%s %lf %lf %d %d %d", temp_algo, &td1, &td2, &t_n, &t_n2, &t_n3);
        }
        /* fprintf(stderr, "%s\n", buff); */
        mcpara->algo=temp_algo;
        mcpara->dfac = td1;
        mcpara->kBT = td2;
        mcpara->is_restart = t_n;
        mcpara->tot_mc_iter = t_n2;
        mcpara->dump_skip = t_n3;
        /* mcpara->activity = td3; */
        if( fgets(buff,255,(FILE*)f2) != NULL);  
        if( fgets(buff,255,(FILE*)f2) != NULL);   
        if( fgets(buff,255,(FILE*)f2) != NULL){   
            sscanf(buff,"%s %lf %lf", which_act, &td1, &td2);
        }
        /* cout << "activity" << "\t"<<f2<<endl; */
       activity->act = which_act;
       activity->minA = td1;
       activity->maxA = td2;
        /* fprintf(stderr, "%s\n", buff); */
       if( fgets(buff,255,(FILE*)f2) != NULL);  
       if( fgets(buff,255,(FILE*)f2) != NULL);   
       if( fgets(buff,255,(FILE*)f2) != NULL){   
           sscanf(buff,"%s %lf %lf", which_act, &td1, &td2);
       }
        afm->tip_rad = td1;
        afm->tip_pos_z = td2;
        afm->sigma = td3;
        afm->epsilon = td4;
        if (fabs(afm->epsilon)<1e-16){
            afm->icompute=0;
        }
        if( fgets(buff,255,(FILE*)f2) != NULL);  
        if( fgets(buff,255,(FILE*)f2) != NULL);   
        if( fgets(buff,255,(FILE*)f2) != NULL){
            sscanf(buff,"%d %lf %lf", &t_n, &td1, &td2);
        }
        spring -> icompute = t_n;
        spring -> nPole_eq_z = td1;
        spring -> sPole_eq_z = td2;

    }

    else{
        fprintf(stderr, "Mayday Mayday the specified para file doesn't exists\n");
        fprintf(stderr, "I will just kill myself now\n");
        exit(0);
    }
    fclose(f2);
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