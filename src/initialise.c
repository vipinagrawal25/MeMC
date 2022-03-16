#include "../include/global.h"
#include "../include/subroutine.h"

int randint(int n) {
    if ((n - 1) == RAND_MAX) {
        return rand();
    } else {
        // Chop off all of the values that would cause skew...
        long end = RAND_MAX / n; // truncate skew
        assert (end > 0L);
        end *= n;
        // ... and ignore results from rand() that fall above that limit.
        // so we can expect to bail out of this loop pretty quickly.)
        int r;
        while ((r = rand()) >= end);

        return r % n;
    }
}

// void initialize_read_config(POSITION *Pos,  
//         MESH mesh,  int *triangles,
//         MBRANE_para mb_para){
//     FILE *fid;
//     int2 *tmp_nbr_lst;

//     tmp_nbr_lst = mesh.bond_nbr_list;
//     // Read the particle postions;
//     //
//     //
//     fid = fopen("analysis/positions.dat", "rb");
//     for(int i=0; i<mb_para.N; i++){
//         fscanf(fid, "%lf %lf %lf \n", &Pos[i].x, &Pos[i].y, &Pos[i].z);
//         /* fprintf(stderr, "%lf %lf %lf \n", Pos[i].x, Pos[i].y, Pos[i].z); */
//     }
//     fclose(fid);

//     fid = fopen("analysis/cumlist.dat", "rb");
//     for(int i=0; i<mb_para.N+1; i++){
//         fscanf(fid, "%d \n", &mesh.cmlist[i]);
//     }
//     fclose(fid);

//     fid = fopen("analysis/node_neighbour.dat", "rb");
//     for(int i=0; i<mesh.cmlist[mb_para.N]; i++){
//         fscanf(fid, "%d \n", &mesh.node_nbr_list[i]);
//     }
//     fclose(fid);


//     fid = fopen("analysis/bond_neighbour.dat", "rb");
//     for(int i=0; i<mesh.cmlist[mb_para.N]; i++){
//         fscanf(fid, "%d %d\n", &tmp_nbr_lst[i].i1, &tmp_nbr_lst[i].i2);
//     }
//     fclose(fid);

//     fid = fopen("analysis/triangles.dat", "rb");
//     for(int i=0; i<mesh.cmlist[mb_para.N]; i=i+3){
//         fscanf(fid, "%d %d %d\n", &triangles[i], &triangles[i+1], &triangles[i+2]);
//     }
//     fclose(fid);
// }

// void initialize_system(POSITION *Pos,  LJpara para, 
//         MCpara mc_para){
//     int root_n = 1;
//     double dr;
//     bool is_sph, is_cart;

//     do{ 
//         root_n++;
//     }while(root_n*root_n < para.N); 

//     /* metric = "random"; */

//     is_sph = false;
//     is_cart = false;
//     if(strcmp(mc_para.metric, "sph") == 0){
//         is_sph = true;
//     }
//     if(strcmp(mc_para.metric, "cart") == 0){
//         is_cart = true;
//     }


//     dr = (double)(para.len)/(double)root_n;
//     dr = dr - dr/64;

//     if(is_cart){
//         for(int i=0; i<para.N; i++){
//             /* if(configuration == "regular"){ */
//                 /* Pos[i].x =  dr*(double)((i%(root_n))%root_n);// + 2*r; */
//                 /* Pos[i].y =  dr*(double)((i/(root_n)));//+ 2*r; */
//             /* } */
//             /* else if(configuration == "random"){ */
//                 Pos[i].x = drand48()*para.len;
//                 Pos[i].y = drand48()*para.len;
//             /* } */
//         }
//         Pos[0].x = drand48()*para.len;
//         Pos[0].y = drand48()*para.len;
//     }

//     if(is_sph){
//         /* printf("Regular arrangement not coded for spherical"); */
//         /* printf("Coordinates \n"); */
//         Pos[0].x = 0;
//         Pos[0].y = 0;
//         Pos[1].x = pi;
//         Pos[1].y = 0;
//         for(int i=2; i<para.N; i++){
//             Pos[i].x = acos(2*drand48() - 1); 
//             Pos[i].y = 2*pi*drand48();
//             /* printf("Coordinates %lf %lf\n", Pos[i].x, Pos[i].y); */
//         }
//         Pos[2].x = acos(2*drand48() - 1); 
//         Pos[2].y = 2*pi*drand48();

//     }
// }

void initialize_eval_lij_t0(POSITION *Pos, MESH mesh, 
        double *lij_t0, MBRANE_para *para){

    double rij;
    POSITION dr;
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

void initialize_afm_tip(AFM_para afm){
    int i, j;
    int idx;
    int iext;
    double dx;
    double x, y;

    iext = (int) sqrt(afm.N);
    dx = (afm.extent[1] - afm.extent[0])/iext;

    for(j=0; j < iext; j++){
        for(i=0; i < iext; i++){
            x = afm.extent[0] + dx*(i + 0.5);
            y = afm.extent[2] + dx*(j + 0.5);
            idx = j*iext + i;
            afm.tip_curve[idx].x = x;
            afm.tip_curve[idx].y = y;
            afm.tip_curve[idx].z = (x/afm.tip_rad)*(x/afm.tip_rad) +
                (y/afm.tip_rad)*(y/afm.tip_rad) 
                + afm.tip_pos_z;
        }

    }

}

void initialize_read_parameters( MBRANE_para *mbrane, 
        AFM_para *afm, MCpara *mcpara, string para_file){
    
    char buff[255];
    int t_n, t_n2, t_n3;
    double td1, td2, td3, td4, td5;
    FILE *f2;
    f2 = fopen(para_file.c_str(), "r");
    if(f2){
        fgets(buff,255,(FILE*)f2); 
        fgets(buff,255,(FILE*)f2); 
        fgets(buff,255,(FILE*)f2);
        sscanf(buff,"%d %lf %lf %lf", &t_n, &td1, &td2, &td3);
        /* fprintf(stderr, "%s\n", buff); */
        mbrane->N = t_n;
        mbrane->coef_bend = td1;
        mbrane->coef_str = td2;
        mbrane->coef_vol_expansion = td3;
        fgets(buff,255,(FILE*)f2);
        fgets(buff,255,(FILE*)f2);
        sscanf(buff,"%lf %lf %lf %lf %lf", &td1,&td2,&td3,&td4,&td5);
        /* fprintf(stderr, "%s\n", buff); */
        mbrane->radius = td1;
        mbrane->pos_bot_wall = td2;
        mbrane->sigma = td3;
        mbrane->epsilon = td4;
        mbrane->sp_curv = td5;
        fgets(buff,255,(FILE*)f2);
        fgets(buff,255,(FILE*)f2); 
        fgets(buff,255,(FILE*)f2); 
        sscanf(buff,"%lf %lf %d %d %d", &td1, &td2, &t_n, &t_n2, &t_n3);
        /* fprintf(stderr, "%s\n", buff); */
        mcpara->dfac = td1;
        mcpara->kBT = td2;
        mcpara->is_restart = t_n;
        mcpara->tot_mc_iter = t_n2;
        mcpara->dump_skip = t_n3;
        fgets(buff,255,(FILE*)f2);
        fgets(buff,255,(FILE*)f2);
        fgets(buff,255,(FILE*)f2); 
        sscanf(buff,"%d %lf %lf %lf %lf", &t_n, &td1, &td2, &td3, &td4);
        /* fprintf(stderr, "%s\n", buff); */
        afm->N = t_n; 
        afm->extent[0] = td1; afm->extent[1]= td2;
        afm->extent[2] = td3; afm->extent[3]= td3;

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
        fprintf(stderr, "sorry man the specified doesn't exists in current dir .\n");
    }
    fclose(f2);
    mbrane->num_triangles = 2*mbrane->N - 4;
    mbrane->num_nbr = 3*mbrane->num_triangles;
    // TODO: compute the average bond length by running over all the lengths.
    // mbrane->av_bond_len = sqrt(8*pi/(2*mbrane->N-4));
   // define the monte carlo parameters
    mcpara->one_mc_iter = 2*mbrane->N;
    mcpara->metric = "sph";
    mcpara->delta = sqrt(8*pi/(2*mbrane->N-4));
}