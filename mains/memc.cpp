#include "../includes/global.h"
#include "../includes/subroutine.h"
//
double start_simulation(Vec3d *Pos, MESH_p mesh, double *lij_t0, 
                    MBRANE_p mbrane_para, MC_p mc_para, STICK_p stick_para,
                    VOL_p vol_para, AFM_p afm_para, ACTIVE_p act_para, 
                    SPRING_p spring_para, FLUID_p fld_para, string outfolder){
    double Pole_zcoord;
    if(!mc_para.is_restart){
        hdf5_io_read_pos( (double *)Pos,  outfolder+"/input.h5");
        hdf5_io_read_mesh((int *) mesh.numnbr,
                (int *) mesh.node_nbr_list, outfolder+"/input.h5");
        init_eval_lij_t0(Pos, mesh, lij_t0,  &mbrane_para, &spring_para, fld_para.is_fluid);
        if(stick_para.do_stick)
            identify_attractive_part(Pos, stick_para.is_attractive, stick_para.theta, mbrane_para.N);
        max(&mesh.nPole,&Pole_zcoord,Pos,mbrane_para.N);
        min(&mesh.sPole,&Pole_zcoord,Pos,mbrane_para.N);
    }else{
        hdf5_io_read_pos( (double *)Pos,  outfolder+"/input.h5");
        hdf5_io_read_mesh((int *) mesh.numnbr,
                (int *) mesh.node_nbr_list, outfolder+"/input.h5");
        init_eval_lij_t0(Pos, mesh, lij_t0,  &mbrane_para, &spring_para, fld_para.is_fluid);
        if(stick_para.do_stick)
            identify_attractive_part(Pos, stick_para.is_attractive, stick_para.theta, mbrane_para.N);
        max(&mesh.nPole,&Pole_zcoord,Pos,mbrane_para.N);
        min(&mesh.sPole,&Pole_zcoord,Pos,mbrane_para.N);
        identify_attractive_part(Pos, stick_para.is_attractive, stick_para.theta, mbrane_para.N);
        hdf5_io_read_pos( (double *)Pos,  outfolder+"/restart.h5");
    }
    return Pole_zcoord;
}

double diag_energies(double *Et, Vec3d *Pos, MESH_p mesh, double *lij_t0, 
        MBRANE_p mbrane_para, MC_p mc_para, STICK_p stick_para,
        VOL_p vol_para, AFM_p afm_para, ACTIVE_p act_para, 
        SPRING_p spring_para, FILE *fid ){
    double vol_sph;
    double Ener_t;
    Vec3d afm_force,spring_force[2];

    /*-----------------------------------------------*/
    /*****  initialize energy values *****/
    Et[0] = stretch_energy_total(Pos, mesh, lij_t0, mbrane_para);
    Et[1] = bending_energy_total(Pos, mesh, mbrane_para);
    fprintf(fid, " %g %g %g", mbrane_para.tot_energy[0], Et[0], Et[1]);
    if(stick_para.do_stick){
        Et[2] = lj_bottom_surf_total(Pos, mbrane_para, stick_para);
        fprintf(fid, " %g", Et[2]);
    }
    if(afm_para.do_afm){
        Et[3] = lj_afm_total(Pos, &afm_force, mbrane_para, afm_para);
        fprintf(fid, " %g", Et[3]);
    }
    if(vol_para.do_volume){
        vol_sph = volume_total(Pos, mesh, mbrane_para);
        double  ini_vol = (4./3.)*pi*pow(mbrane_para.radius,3);
        Et[4] = vol_para.coef_vol_expansion*(vol_sph/ini_vol - 1e0)*(vol_sph/ini_vol - 1e0);
        mbrane_para.volume[0] = vol_sph;
        fprintf(fid, " %g", Et[4]);
        if(vol_para.is_pressurized){
            Et[6] = -vol_para.pressure*vol_sph;
            fprintf(fid, " %g", Et[6]);
        }
    }
    if(spring_para.do_spring){
        Et[5] = spring_tot_energy_force(Pos, spring_force, mesh, spring_para);
        fprintf(fid, " %g\n", Et[5]);
    }
    Ener_t = Et[0] + Et[1] + Et[2] + Et[3] + Et[4]+ Et[5]+ Et[6];
    fflush(fid);

    return Ener_t;
}

int main(int argc, char *argv[]){
    pid_t pid = getpid();
    cout << "# ID for this process is: " << pid << endl;
    //
    int i, iterations, num_moves, num_bond_change;
    double Et[7], Ener_t;
    double vol_sph, e_t, s_t;
    Vec3d *Pos;
    MBRANE_p mbrane_para;
    MC_p mc_para;
    AFM_p  afm_para;
    ACTIVE_p act_para;
    MESH_p mesh;
    VOL_p vol_para;
    STICK_p stick_para;
    SPRING_p spring_para;
    FLUID_p fld_para;
    Vec3d afm_force,spring_force[2];
    FILE *fid;
    double *lij_t0;
    double Pole_zcoord;
    string outfolder,syscmds, para_file, log_file, outfile, filename;
    int mpi_err,mpi_rank=0;
    int iter;
    uint32_t seed_v;

    //
    mpi_err = MPI_Init(0x0, 0x0);
    mpi_err =  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    seed_v = (uint32_t) 7*3*11*(mpi_rank+1)*rand();
    init_rng(seed_v);
    //
    outfolder = ZeroPadNumber(mpi_rank)+"/";
    cout << "I am in folder "+ outfolder << endl;
    filename = outfolder + "/para_file.in";

    // ---------- open outfile_terminal ------------------- //
    fstream outfile_terminal(outfolder+"/terminal.out", ios::app);
    /*************************************************/
    // read the input file
    init_read_parameters(&mbrane_para, &mc_para, &fld_para, &vol_para,
            &stick_para, &afm_para,  &act_para, &spring_para, filename);

    mc_para.one_mc_iter = 2*mbrane_para.N;


   // check whether the string comparison works
   /* define all the paras */
    mbrane_para.volume = (double *)calloc(1, sizeof(double)); 
    mbrane_para.volume[0] = (4./3.)*pi*pow(mbrane_para.radius,3);
    mbrane_para.tot_energy = (double *)calloc(1, sizeof(double));
    act_para.activity = (double *)calloc(mbrane_para.N, sizeof(double));
    mbrane_para.tot_energy[0] = 0e0;
    init_activity(act_para, mbrane_para.N);
    // allocate arrays
    Pos = (Vec3d *)calloc(mbrane_para.N, sizeof(Vec3d));
    mesh.numnbr = (int *)calloc(mbrane_para.N, sizeof(int));
    mesh.nghst = 12;
    mbrane_para.len = 1.e0;
    mesh.node_nbr_list = (int *)calloc(mesh.nghst*mbrane_para.N, sizeof(int));
    lij_t0 = (double *)calloc(mesh.nghst*mbrane_para.N, sizeof(double));
    stick_para.is_attractive = (bool *)calloc(mbrane_para.N, sizeof(bool));

    //
    if(!mc_para.is_restart && afm_para.do_afm){
        s_t = afm_para.sigma; 
        afm_para.sigma = 0.00;
        e_t = afm_para.epsilon;
        afm_para.epsilon = 0.0;
    }

    Pole_zcoord = start_simulation(Pos, mesh, lij_t0, 
                     mbrane_para,  mc_para,  stick_para,
                     vol_para,  afm_para,  act_para, 
                     spring_para,  fld_para,  outfolder);
    /************************************/
    // cout << "# Foppl von Karman (FvK): " 
         // << YY*mbrane.radius*mbrane.radius/BB << endl;
    //
    log_file=outfolder+"/mc_log";
    fid = fopen(log_file.c_str(), "a");
    /* wHeader(fid,mbrane,afm,spring); */
    Ener_t = diag_energies(Et, Pos,  mesh, lij_t0,  mbrane_para,  mc_para,  stick_para,
         vol_para,  afm_para,  act_para, spring_para,  fid );
    num_moves = 0;
    for(i=0; i < mc_para.tot_mc_iter; i++){

        if(i%mc_para.dump_skip == 0){
            outfile=outfolder+"/snap_"+ZeroPadNumber(i/mc_para.dump_skip)+".h5";
            hdf5_io_write_pos((double*) Pos, 3*mbrane_para.N, outfile);
            syscmds="cp "+outfile+" "+outfolder+"/restart.h5";
            system(syscmds.c_str());
        }
        if(i == 10*mc_para.dump_skip && !mc_para.is_restart && afm_para.do_afm){
            afm_para.sigma = s_t;
            afm_para.epsilon = e_t;
            e_t = lj_afm_total(Pos, &afm_force, mbrane_para, afm_para);
            mbrane_para.tot_energy[0] += e_t;
        }

        num_moves = monte_carlo_3d(Pos, mesh, lij_t0, 
                mbrane_para, mc_para, stick_para, vol_para, 
                afm_para, act_para,  spring_para);

        if(fld_para.is_fluid && iter%fld_para.fluidize_every==0){
            num_bond_change = monte_carlo_fluid(Pos, mesh, mbrane_para, mc_para, fld_para);
        }

        fprintf(fid , "%d %d", i, num_moves );
        Ener_t = diag_energies(Et, Pos,  mesh, lij_t0,  mbrane_para,  mc_para,  stick_para,
                vol_para,  afm_para,  act_para, spring_para,  fid );
        outfile_terminal << "iter = " << i << "; Accepted Moves = " 
            << (double) num_moves*100/mc_para.one_mc_iter << " %;"<<  
            " totalener = "<< mbrane_para.tot_energy[0] << "; volume = " << vol_sph << endl;

    }

    fclose(fid);
    free(Pos);
    free(lij_t0);
    free(mesh.node_nbr_list);
    mpi_err = MPI_Finalize();
    outfile_terminal.close();
    return 0;
}
