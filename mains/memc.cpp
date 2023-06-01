#include "../includes/global.h"
#include "../includes/subroutine.h"
//
int main(int argc, char *argv[]){
    pid_t pid = getpid();
    cout << "# ID for this process is: " << pid << endl;
    //
    int i, iterations, num_moves;
    double Et[7], Ener_t;
    double vol_sph, e_t, s_t;
    Vec3d *Pos;
    bool *is_attractive;
    MBRANE_para mbrane;
    MCpara mcpara;
    AFM_para afm;
    ActivePara activity;
    MESH mesh, mes_t;
    Vec3d afm_force,spring_force[2];
    FILE *fid;
    double *lij_t0;
    double Pole_zcoord;
    string outfolder,syscmds, para_file, log_file, outfile, filename;
    int mpi_err,mpi_rank;
    uint32_t seed_v;
    char log_headers[] = "# iter acceptedmoves total_ener stretch_ener bend_ener stick_ener afm_ener ener_volume";
    SPRING_para spring;
    //
    mpi_err = MPI_Init(0x0, 0x0);
    mpi_err =  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    seed_v = (uint32_t) 7*3*11*(mpi_rank+1)*rand();
    init_rng(seed_v);
    //
    outfolder = ZeroPadNumber(mpi_rank)+"/";
    // cout << "I am in folder "+ outfolder << endl;
    filename = outfolder + "/para_file.in";
    // ---------- open outfile_terminal ------------------- //
    fstream outfile_terminal(outfolder+"/terminal.out", ios::app);
    /*************************************************/
    // read the input file
    init_read_parameters(&mbrane, &afm, &mcpara, &activity, &spring, filename);
    /* define all the paras */ 
    mbrane.volume = (double *)calloc(1, sizeof(double)); 
    mbrane.volume[0] = (4./3.)*pi*pow(mbrane.radius,3);
    mbrane.tot_energy = (double *)calloc(1, sizeof(double));
    activity.activity = (double *)calloc(mbrane.N, sizeof(double));
    mbrane.tot_energy[0] = 0e0;
    init_activity(activity, mbrane.N);
    write_param(outfolder + "/para.out",mbrane,mcpara,spring);
    // allocate arrays
    Pos = (Vec3d *)calloc(mbrane.N, sizeof(Vec3d));
    mesh.cmlist = (int *)calloc(mbrane.N+1, sizeof(int));
    mesh.node_nbr_list = (int *)calloc(mbrane.num_nbr, sizeof(int)); 
    lij_t0 = (double *)calloc(mbrane.num_nbr, sizeof(double));
    is_attractive = (bool *)calloc(mbrane.N, sizeof(bool));
    //
    if(!mcpara.is_restart){
        s_t = afm.sigma; 
        afm.sigma = 0.00;
        e_t = afm.epsilon;
        afm.epsilon = 0.0;
    }
    //
    if(!mcpara.is_restart){
        hdf5_io_read_pos( (double *)Pos,  outfolder+"/dmemc_pos.h5");
        hdf5_io_read_mesh((int *) mesh.cmlist,
                (int *) mesh.node_nbr_list, outfolder+"/dmemc_conf.h5");
        init_eval_lij_t0(Pos, mesh, lij_t0, &mbrane, &spring);
        identify_attractive_part(Pos, is_attractive, mbrane.theta, mbrane.N);
        max(&mesh.nPole,&Pole_zcoord,Pos,mbrane.N);
        min(&mesh.sPole,&Pole_zcoord,Pos,mbrane.N);
    }else{
        hdf5_io_read_pos( (double *)Pos, outfolder+"/dmemc_pos.h5");
        hdf5_io_read_mesh((int *) mesh.cmlist,
                (int *) mesh.node_nbr_list,  outfolder+"/dmemc_conf.h5");
        max(&mesh.nPole,&Pole_zcoord,Pos,mbrane.N);
        min(&mesh.sPole,&Pole_zcoord,Pos,mbrane.N);
        init_eval_lij_t0(Pos, mesh, lij_t0, &mbrane, &spring);
        identify_attractive_part(Pos, is_attractive, mbrane.theta, mbrane.N);
        hdf5_io_read_pos( (double *)Pos,  outfolder+"/restart.h5");
    }
    // int j,k;
    // int num_nbr, cm_idx;
    // for(i = 0; i < mbrane.N; i++){
    //     num_nbr = mesh.cmlist[i + 1] - mesh.cmlist[i];
    //     cm_idx = mesh.cmlist[i];
    //     for(k = cm_idx; k < cm_idx + num_nbr; k++) {
    //         cout << lij_t0[k] << "\t";
    //     }
    // }
    // cout << endl;
    /*****  initialize energy values *****/
    Et[0] = stretch_energy_total(Pos, mesh, lij_t0, mbrane);
    Et[1] = bending_energy_total(Pos, mesh, mbrane);
    Et[2] = lj_bottom_surf_total(Pos, is_attractive, mbrane);
    Et[3] = lj_afm_total(Pos, &afm_force, mbrane, afm);
    vol_sph = volume_total(Pos, mesh, mbrane);
    double  ini_vol = (4./3.)*pi*pow(mbrane.radius,3);
    Et[4] = mbrane.coef_vol_expansion*(vol_sph/ini_vol - 1e0)*(vol_sph/ini_vol - 1e0);
    Et[5] = spring_tot_energy_force(Pos, spring_force, mesh, spring);
    Et[6] = -mbrane.pressure*vol_sph;
    Ener_t = Et[0] + Et[1] + Et[2] + Et[3] + Et[4]+ Et[5]+ Et[6];
    mbrane.tot_energy[0] = Ener_t;
    mbrane.volume[0] = vol_sph;
    // cout << Et[0] << "\t" << Et[1] << "\t" << vol_sph  << "\t" << Ener_t << endl;
    /************************************/
    log_file=outfolder+"/mc_log";
    fid = fopen(log_file.c_str(), "a");
    wHeader(fid,mbrane,afm,spring);
    num_moves = 0;
    for(i=0; i < mcpara.tot_mc_iter; i++){
        Et[0] =  stretch_energy_total(Pos, mesh, lij_t0, mbrane);
        Et[1] =  bending_energy_total(Pos, mesh, mbrane);
        Et[2] = lj_bottom_surf_total(Pos, is_attractive, mbrane);
        Et[3] = lj_afm_total(Pos, &afm_force, mbrane, afm);
        vol_sph = volume_total(Pos, mesh, mbrane);
        Et[4] = mbrane.coef_vol_expansion*(vol_sph/ini_vol - 1e0)*(vol_sph/ini_vol - 1e0);
        Et[5] = spring_tot_energy_force(Pos, spring_force, mesh, spring);
        Et[6] = -mbrane.pressure*vol_sph;
        // cout << "iter = " << i << "; Accepted Moves = " << (double) num_moves*100/mcpara.one_mc_iter << " %;"<<  
        //         " totalener = "<< mbrane.tot_energy[0] << "; volume = " << vol_sph << endl;
        outfile_terminal << "iter = " << i << "; Accepted Moves = " << (double) num_moves*100/mcpara.one_mc_iter << " %;"<<  
                " totalener = "<< mbrane.tot_energy[0] << "; volume = " << vol_sph << endl;
        wDiag(fid, mbrane, afm, spring, mesh, i, num_moves, Et,  &afm_force,  spring_force,
              vol_sph, Pos);
        if(i%mcpara.dump_skip == 0){
            outfile=outfolder+"/snap_"+ZeroPadNumber(i/mcpara.dump_skip)+".h5";
            // sprintf(outfile,"%s/snap_%04d.h5",outfolder,(int)(i/mcpara.dump_skip));
            hdf5_io_write_pos((double*) Pos, 3*mbrane.N, outfile);
            syscmds="cp "+outfile+" "+outfolder+"/restart.h5";
            system(syscmds.c_str());
        }
        if(i == 10*mcpara.dump_skip && !mcpara.is_restart){
            afm.sigma = s_t;
            afm.epsilon = e_t;
            e_t = lj_afm_total(Pos, &afm_force, mbrane, afm);
            mbrane.tot_energy[0] += e_t;
        }
        num_moves = monte_carlo_3d(Pos, mesh, lij_t0, is_attractive, 
                mbrane, mcpara, afm, activity,  spring);
    }
    fclose(fid);
    free(Pos);
    free(lij_t0);
    free(mesh.node_nbr_list);
    mpi_err = MPI_Finalize();
    outfile_terminal.close();
    return 0;
}