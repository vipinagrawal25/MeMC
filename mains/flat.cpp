#include "../includes/global.h"
#include "../includes/subroutine.h"
#include<assert.h>

void check_absurd_neighbours(MESH mesh, int N){
    int i, cmidx, nnbr, j;
    for(i=0;i < N; i++){
            cmidx = i*mesh.nghst;
            nnbr = mesh.numnbr[i];
            for(j=0; j< nnbr; j++ ){
                assert(mesh.node_nbr_list[cmidx+j] != -1);
            }

        }
}

void debug_ener(Vec3d *pos, MESH mesh, 
                double *lij_t0,  MBRANE_para mbrane, 
                MCpara mcpara, AFM_para afm, 
                ActivePara activity, SPRING_para spring){

    int idx;
    int cm_idx, num_nbr;
    double E_b;

    for(idx = 0; idx < mbrane.N; idx++){
        cm_idx = mesh.nghst*idx;
        num_nbr = mesh.numnbr[idx];


        /* E_b = bending_energy_ipart_neighbour(pos, mesh, idx, mbrane); */
      E_b = stretch_energy_ipart(pos, (int *) (mesh.node_nbr_list + cm_idx),
                 (double *) (lij_t0 + cm_idx), num_nbr,
                 idx, mbrane);


        /* E_b = bending_energy_ipart(pos, */ 
        /*         (int *) (mesh.node_nbr_list + cm_idx), */
        /*         num_nbr, idx, mbrane); */

        printf(" curvature %d %lf \n",idx, E_b);
    }
}
//
int main(int argc, char *argv[]){
    pid_t pid = getpid();
    cout << "# ID for this process is: " << pid << endl;
    //
    int i, iterations, num_moves;
    double Et[7], Ener_t;
    Vec3d *Pos;
    MBRANE_para mbrane;
    MCpara mcpara;
    AFM_para afm;
    ActivePara activity;
    MESH mesh;
    Vec3d afm_force,spring_force[2];
    FILE *fid;
    double *lij_t0;
    bool *is_attractive;
    double Pole_zcoord;
    string outfolder,syscmds, para_file, log_file, outfile, filename;
    int mpi_err,mpi_rank=0;
    int iter;
    uint32_t seed_v;
    char log_headers[] = "# iter acceptedmoves total_ener stretch_ener bend_ener stick_ener afm_ener ener_volume";
    SPRING_para spring;
    //
    /* mpi_err = MPI_Init(0x0, 0x0); */
    /* mpi_err =  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); */
    /* seed_v = (uint32_t) 7*3*11*(mpi_rank+1)*rand(); */
    seed_v = 12343234;
    init_rng(seed_v);
    //
    outfolder = argv[1]; //Z"e"roPadNumber(mpi_rank)+"/";
    cout << "I am in folder "+ outfolder << endl;
    filename = outfolder + "/para_file.in";
    write_param(outfolder + "/para.out",mbrane,mcpara,spring);
        // ---------- open outfile_terminal ------------------- //
    fstream outfile_terminal(outfolder+"/terminal.out", ios::app);
    /*************************************************/
    // read the input file
    init_read_parameters(&mbrane, &afm, &mcpara, &activity, 
                        &spring, filename);

    mcpara.one_mc_iter = 2*mbrane.N;
    // keep incase I want to check some parameters
   fprintf(stderr, "%d, %g %g %g \n", mbrane.N, mbrane.YY, mbrane.coef_bend, mbrane.coef_vol_expansion );
   fprintf(stderr, "%g, %g %g %g \n", mbrane.pressure, mbrane.sigma, mbrane.epsilon, mbrane.theta );
   fprintf(stderr, "%d, %g %g %g \n", mbrane.N, mbrane.YY, mbrane.coef_bend, mbrane.coef_vol_expansion );
   fprintf(stderr, "%s, %d %g %g \n", activity.act.c_str(), mcpara.one_mc_iter, activity.maxA, mbrane.coef_vol_expansion );

   // check whether the string comparison works
   /* define all the paras */
    mbrane.tot_energy = (double *)calloc(1, sizeof(double));
    mbrane.len = 1.0;
    mbrane.is_fluid = true;
    activity.activity = (double *)calloc(mbrane.N, sizeof(double));
    mbrane.tot_energy[0] = 0e0;
    init_activity(activity, mbrane.N);
    // allocate arrays
    Pos = (Vec3d *)calloc(mbrane.N, sizeof(Vec3d));
    mesh.numnbr = (int *)calloc(mbrane.N, sizeof(int));
    mesh.nghst = 12;
    mesh.node_nbr_list = (int *)calloc(mesh.nghst*mbrane.N, sizeof(int));
   lij_t0 = (double *)calloc(mesh.nghst*mbrane.N, sizeof(double));


    // if(!mcpara.is_restart){
        hdf5_io_read_pos( (double *)Pos,  outfolder+"/inp_pos.h5");
        hdf5_io_read_mesh((int *) mesh.numnbr,
                (int *) mesh.node_nbr_list, outfolder+"/nbrs.h5");
        init_eval_lij_t0(Pos, mesh, lij_t0,  &mbrane, &spring);
    // Et[0] = stretch_energy_total(Pos, mesh, lij_t0, mbrane);
    


    for(iter=1; iter <= mcpara.tot_mc_iter; iter++){

     num_moves = monte_carlo_3d(Pos, mesh, lij_t0, is_attractive,
           mbrane, mcpara, afm, activity,  spring);
          
    printf("elastic part stats %d %d \n", iter, num_moves );
     /* num_moves = monte_carlo_fluid(Pos, mesh, */
     /*         mbrane, mcpara, afm, activity,  spring); */
    /* printf("fluid stats %d %d \n", iter, num_moves ); */
     if(i%mcpara.dump_skip == 0){
         outfile=outfolder+"/snap_"+ZeroPadNumber(iter/mcpara.dump_skip)+".h5";
         // sprintf(outfile,"%s/snap_%04d.h5",outfolder,(int)(i/mcpara.dump_skip));
         hdf5_io_write_pos((double*) Pos, 3*mbrane.N, outfile);
         /* syscmds="cp "+outfile+" "+outfolder+"/restart.h5"; */
         system(syscmds.c_str());
     }
 
    /* check_absurd_neighbours(mesh, mbrane.N); */

    }

    // /************************************/
    // // cout << "# Foppl von Karman (FvK): "
    //      // << YY*mbrane.radius*mbrane.radius/BB << endl;
    // //

    /* fid = fopen(log_file.c_str(), "a"); */

    //     Et[0] =  stretch_energy_total(Pos, mesh, lij_t0, mbrane);
    //     cout << "iter = " << i << "; Accepted Moves = " << (double) num_moves*100/mcpara.one_mc_iter << " %;"<<
    //             " totalener = "<< mbrane.tot_energy[0] << "; volume = " << vol_sph << endl;
   //     num_moves = monte_carlo_fluid(Pos, mesh, mbrane, mcpara, afm, activity,  spring);

    // }
    /* fclose(fid); */
    free(Pos);
    free(mesh.node_nbr_list);
    /* mpi_err = MPI_Finalize(); */
    return 0;
}
