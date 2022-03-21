#include "global.h"
#include "subroutine.h"

void init_system_random_pos(Vec2d *Pos,  double len, int N, char *metric){
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
        double *lij_t0, MBRANE_para *para){

    Vec3d dr;
    int i,j,k;
    int num_nbr, cm_idx;
    double sum_lij=0;
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
}

void init_read_parameters( MBRANE_para *mbrane, 
        AFM_para *afm, MCpara *mcpara, char *para_file){
    
    char buff[255];
    int t_n, t_n2, t_n3;
    double td1, td2, td3, td4, td5;
    FILE *f2;
    f2 = fopen(para_file, "r");
    if(f2){
        fgets(buff,255,(FILE*)f2); 
        fgets(buff,255,(FILE*)f2); 
        fgets(buff,255,(FILE*)f2);
        sscanf(buff,"%d %lf %lf %lf", &t_n, &td1, &td2, &td3, &td4);
        /* fprintf(stderr, "%s\n", buff); */
        mbrane->N = t_n;
        mbrane->coef_bend = td1;
        mbrane->coef_str = td2;
        mbrane->coef_vol_expansion = td3;
        mbrane->sp_curv = td4;
        fgets(buff,255,(FILE*)f2);
        fgets(buff,255,(FILE*)f2);
        sscanf(buff,"%lf %lf %lf %lf %lf", &td1,&td2,&td3,&td4,&td5);
        /* fprintf(stderr, "%s\n", buff); */
        mbrane->radius = td1;
        mbrane->pos_bot_wall = td2;
        mbrane->sigma = td3;
        mbrane->epsilon = td4;
        mbrane->theta = td5;
        fgets(buff,255,(FILE*)f2);
        fgets(buff,255,(FILE*)f2); 
        fgets(buff,255,(FILE*)f2); 
        sscanf(buff,"%lf %lf %d %d %d", &td1, &td2,  &t_n2, &t_n3);
        /* fprintf(stderr, "%s\n", buff); */
        mcpara->dfac = td1;
        mcpara->kBT = td2;
        mcpara->tot_mc_iter = t_n2;
        mcpara->dump_skip = t_n3;
        fgets(buff,255,(FILE*)f2);
        fgets(buff,255,(FILE*)f2); 
        fgets(buff,255,(FILE*)f2); 
        sscanf(buff,"%lf %lf %lf %lf", &td1, &td2, &td3, &td4);
        /* fprintf(stderr, "%s\n", buff); */
        afm->tip_rad = td1;
        afm->tip_pos_z = td2;
        afm->sigma = td3;
        afm->epsilon = td4;
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
