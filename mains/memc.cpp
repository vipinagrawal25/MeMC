#include "Vector.h"
#include "global.h"
#include "subroutine.h"
#include <cstring>
//
int * start_simulation(Vec3d *Pos, Vec3d *Pos_t0, MESH_p mesh, double *lij_t0, 
                    MBRANE_p mbrane_para, MC_p mc_para, STICK_p stick_para,
                    VOL_p vol_para, AFM_p afm_para, ACTIVE_p act_para, 
                     FLUID_p fld_para, string outfolder){
    double Pole_zcoord;
    int* poleidx = new int[2];
    if(!mc_para.is_restart){
        hdf5_io_read_double( (double *)Pos,  outfolder+"/input.h5", "pos");
        memcpy((double*) Pos_t0, (double *) Pos, sizeof(double)* 3*mbrane_para.N);
        hdf5_io_read_mesh((int *) mesh.numnbr,
                (int *) mesh.node_nbr_list, outfolder+"/input.h5");
        init_eval_lij_t0(Pos_t0, mesh, lij_t0,  &mbrane_para,  fld_para.is_fluid);
        init_stick_bottom_new(Pos, mesh, stick_para, fld_para, mbrane_para, outfolder );
        hdf5_io_dump_int(fld_para.solid_idx, mbrane_para.N, outfolder+ "/solid_pos_id.h5", "solid_id");
        hdf5_io_dump_bool(stick_para.is_attractive, mbrane_para.N, outfolder+ "/stick_pos_id.h5", "stick_id");
        max(&mesh.nPole,&Pole_zcoord,Pos,mbrane_para.N);
        min(&mesh.sPole,&Pole_zcoord,Pos,mbrane_para.N);
        poleidx[0] = mesh.nPole;
        poleidx[1] = mesh.sPole;
        /* exit(0); */
    }else{
        hdf5_io_read_double( (double *)Pos,  outfolder+"/input.h5", "pos");
        memcpy((double*) Pos_t0, (double *) Pos, sizeof(double)* 3*mbrane_para.N);
        hdf5_io_read_mesh((int *) mesh.numnbr,
                (int *) mesh.node_nbr_list, outfolder+"/input.h5");
        init_eval_lij_t0(Pos, mesh, lij_t0,  &mbrane_para, fld_para.is_fluid);
        // if(stick_para.do_stick)
            // identify_attractive_part(Pos, stick_para.is_attractive, stick_para.theta, mbrane_para.N);
        max(&mesh.nPole,&Pole_zcoord,Pos,mbrane_para.N);
        min(&mesh.sPole,&Pole_zcoord,Pos,mbrane_para.N);
        // if(stick_para.do_stick)
            // identify_attractive_part(Pos, stick_para.is_attractive, stick_para.theta, mbrane_para.N);
        hdf5_io_read_double( (double *)Pos,  outfolder+"/restart.h5", "pos");
        hdf5_io_read_mesh((int *) mesh.numnbr,
                (int *) mesh.node_nbr_list, outfolder+"/restart.h5");
    }
    return poleidx;
}

void diag_wHeader(MBRANE_p mbrane_para, AREA_p area_para, STICK_p stick_para,
        VOL_p vol_para, AFM_p afm_para, ACTIVE_p act_para, 
        SHEAR_p shear_para, FILE *fid ){

    string log_headers = "#iter acceptedmoves total_e bend_e stretch_e ";
    /* if(area_para.is_stick){log_headers+="stick_e ";} */
    if(stick_para.do_stick){log_headers+="stick_e ";}
    if(afm_para.do_afm){log_headers+="afm_e ";}
    if (shear_para.do_shear){log_headers+="spring_e ";}
    if(vol_para.do_volume && !vol_para.is_pressurized) {log_headers+="vol_e ";}
    if (vol_para.do_volume && vol_para.is_pressurized){log_headers+="pressure_e ";}
    if (afm_para.do_afm){log_headers+="afm_fx, afm_fy afm_fz ";}
    log_headers+=" volume  area";
    fprintf(fid, "%s\n", log_headers.c_str());
    fflush(fid);
}

double diag_energies(double *Et, Vec3d *Pos, Vec3d *Pos_t0, MESH_p mesh, double *lij_t0, 
        MBRANE_p mbrane_para, AREA_p area_para, STICK_p stick_para,
        VOL_p vol_para, AFM_p afm_para, ACTIVE_p act_para, 
         FILE *fid ){
    double vol_sph, ar_sph, ini_ar;
    double Ener_t;
    Vec3d afm_force,spring_force[2];

    /*-----------------------------------------------*/
    /*****  initialize energy values *****/

   Et[0] = 0e0;    Et[1] = 0e0;     Et[2] = 0e0; Et[3] = 0e0; 
   Et[4] = 0e0;    Et[5] = 0e0;     Et[6] = 0e0; 

    Et[0] = bending_energy_total(Pos, mesh, mbrane_para);
    ar_sph = area_total(Pos, mesh, mbrane_para);
    /* *mbrane_para.area = ar_sph; */
    if(area_para.is_tether){
        Et[1] = stretch_energy_total(Pos, mesh, lij_t0, mbrane_para, area_para);
    }else{
        ini_ar = (4.)*pi*pow(mbrane_para.radius,2);
        Et[1] = area_para.Ka*(ar_sph/ini_ar - 1e0)*(ar_sph/ini_ar - 1e0);
    }
    fprintf(fid, " %g %g %g", *mbrane_para.tot_energy, Et[0], Et[1]);
    if(stick_para.do_stick){
        Et[2] = lj_bottom_surf_total(Pos, Pos_t0, mbrane_para, stick_para);
        fprintf(fid, " %g", Et[2]);
    }
    if(afm_para.do_afm){
        Et[3] = lj_afm_total(Pos, &afm_force, mbrane_para, afm_para);
        fprintf(fid, " %g", Et[3]);
    }
    
    vol_sph = volume_total(Pos, mesh, mbrane_para);
    /* *mbrane_para.volume = vol_sph; */
    if(vol_para.do_volume){
        double  ini_vol = (4./3.)*pi*pow(mbrane_para.radius,3);
        if(!vol_para.is_pressurized){
            Et[4] = vol_para.coef_vol_expansion*(vol_sph/ini_vol - 1e0)*(vol_sph/ini_vol - 1e0);
            fprintf(fid, " %g", Et[4]);
        }
        if(vol_para.is_pressurized){
            Et[6] = -vol_para.pressure*vol_sph;
            fprintf(fid, " %g", Et[6]);
        }
    }
    if (afm_para.do_afm){fprintf(fid, " %g %g %g", afm_force.x, afm_force.y, afm_force.z);}
    {fprintf(fid, " %g  %g\n", vol_sph, ar_sph );}
    Ener_t = Et[0] + Et[1] + Et[2] + Et[3] + Et[4]+ Et[5]+ Et[6];
    fflush(fid);

    return Ener_t;
}

int main(int argc, char *argv[]){
    pid_t pid = getpid();
    int iter, num_moves, num_bond_change;
    double Et[7], Ener_t;
    double vol_sph, e_t, s_t;
    Vec3d *Pos; MBRANE_p mbrane_para;
    Vec3d *Pos_t0; 
    bool *mask_ids;
    MC_p mc_para; AFM_p  afm_para;
    ACTIVE_p act_para; MESH_p mesh;
    VOL_p vol_para; STICK_p stick_para;
    SHEAR_p shear_para; FLUID_p fld_para;
    AREA_p area_para;
    Vec3d afm_force,spring_force[2];
    Vec2d eners;
    FILE *fid, *fp2;
    double *lij_t0;
    double Pole_zcoord;
    double start_time, end_time;
    string outfolder,syscmds, ener_file, log_file, outfile, filename;
    int mpi_err,mpi_rank=0;
    uint32_t seed_v, num_shear_moves;


    //
    mpi_err = MPI_Init(0x0, 0x0);
    mpi_err =  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    /* seed_v = 12345; //(uint32_t) 7*3*11*(mpi_rank+1)*rand(); */
    seed_v = (uint32_t) (mpi_rank + time(0));

    /* init_rng(seed_v); */
    /* init_rng2(seed_v); */
    RandomGenerator::init(seed_v);
    //
    //

    cout << "# ID for this process is: " << pid << endl;
    outfolder = ZeroPadNumber(mpi_rank+atoi(argv[1]))+"/";
    cout << "# I am in folder "+ outfolder << endl;
    filename = outfolder + "/para_file.in";

    // ---------- open outfile_terminal ------------------- //
    fstream outfile_terminal(outfolder+"/terminal.out", ios::app);
    outfile_terminal << "The seed is  " << seed_v << endl; 
    /*************************************************/

    // read the input file
    init_read_parameters(&mbrane_para, &mc_para, &area_para, &fld_para, &vol_para,
            &stick_para, &afm_para,  &act_para, &shear_para, filename);

    mc_para.one_mc_iter = 2*mbrane_para.N;


   // check whether the string comparison works
   /* define all the paras */
    mbrane_para.volume =  (double*)malloc(sizeof(double));
    *mbrane_para.volume = (4./3.)*pi*pow(mbrane_para.radius,3);
    act_para.activity = (double *)calloc(mbrane_para.N, sizeof(double));
    mbrane_para.tot_energy =  (double*)malloc(sizeof(double));
    mbrane_para.area =  (double*)malloc(sizeof(double));
    *mbrane_para.tot_energy = 0e0;
    init_activity(act_para, mbrane_para.N);
    // allocate arrays
    Pos = (Vec3d *)calloc(mbrane_para.N, sizeof(Vec3d));
    Pos_t0 = (Vec3d *)calloc(mbrane_para.N, sizeof(Vec3d));
    mesh.numnbr = (int *)calloc(mbrane_para.N, sizeof(int));
    mesh.nghst = 12;
    mbrane_para.len = 2.e0*pi;
    mesh.node_nbr_list = (int *)calloc(mesh.nghst*mbrane_para.N, sizeof(int));
    lij_t0 = (double *)calloc(mesh.nghst*mbrane_para.N, sizeof(double));
    stick_para.is_attractive = (bool *)calloc(mbrane_para.N, sizeof(bool));
    fld_para.solid_idx = (int *)calloc(mbrane_para.N, sizeof(int));
    mask_ids = (bool *)calloc(mbrane_para.N, sizeof(bool));

    //
    if(!mc_para.is_restart && afm_para.do_afm){
        s_t = afm_para.sigma; 
        afm_para.sigma = 0.00;
        e_t = afm_para.epsilon;
        afm_para.epsilon = 0.0;
    }
   int *poleidx;

   poleidx = start_simulation(Pos, Pos_t0, mesh, lij_t0, 
                     mbrane_para,  mc_para,  stick_para,
                     vol_para,  afm_para,  act_para, 
                     fld_para,  outfolder);
   mesh.nPole = poleidx[0];
   mesh.sPole = poleidx[1];
    //
    //

    mask_frame(mask_ids, mesh, mbrane_para);
    if(fld_para.is_fluid)mbrane_para.av_bond_len = lij_t0[0];
    log_file=outfolder+"/mc_log";
    ener_file=outfolder+"/ener_log";
    fid = fopen(log_file.c_str(), "a");
    fp2 = fopen(ener_file.c_str(), "a");

    if(!mc_para.is_restart)
    diag_wHeader(mbrane_para, area_para,  stick_para,  vol_para,  afm_para,  act_para, 
         shear_para, fid );



    fprintf(fid , "%d %g", 0, 0.0 );
    Ener_t = diag_energies(Et, Pos, Pos_t0, mesh, lij_t0,  mbrane_para, area_para,  stick_para,
         vol_para,  afm_para,  act_para,   fid );
    *mbrane_para.tot_energy = Ener_t;
    filename = outfolder + "/para.out";
    write_parameters(mbrane_para, mc_para, area_para, fld_para, vol_para,
            stick_para, afm_para,  act_para, shear_para, filename);

    double ar_sph = area_total(Pos, mesh, mbrane_para);
    *mbrane_para.area = ar_sph;
    vol_sph = volume_total(Pos, mesh, mbrane_para);
    *mbrane_para.volume = vol_sph;
    //printf("%lf \n", mbrane_para.tot_energy[0]);
    num_moves = 0;
    start_time = MPI_Wtime();
    double shear_slope = 0.0;

    for(iter=0; iter < mc_para.tot_mc_iter; iter++){

        if(iter%mc_para.dump_skip == 0){
            outfile=outfolder+"/snap_"+ZeroPadNumber(iter/mc_para.dump_skip)+".h5";
            hdf5_io_write_double((double*) Pos, 3*mbrane_para.N, outfile, "pos");
            hdf5_io_write_mesh(mesh.numnbr, mesh.node_nbr_list, 
                    mbrane_para.N, mesh.nghst, outfile);
            syscmds="cp "+outfile+" "+outfolder+"/restart.h5";
            system(syscmds.c_str());
        }
        if(iter == 10*mc_para.dump_skip && !mc_para.is_restart && afm_para.do_afm){
            afm_para.sigma = s_t;
            afm_para.epsilon = e_t;
            e_t = lj_afm_total(Pos, &afm_force, mbrane_para, afm_para);
            mbrane_para.tot_energy[0] += e_t;
        }

        num_moves = monte_carlo_3d(Pos, Pos_t0, mesh, lij_t0, 
                mbrane_para, mc_para, area_para, stick_para, vol_para, 
                afm_para, act_para,  shear_para);

        if(fld_para.is_fluid && iter%fld_para.fluidize_every==0){
            num_bond_change = monte_carlo_fluid(Pos, mesh, mbrane_para, mc_para, fld_para);
            outfile_terminal << "fluid stats " << num_bond_change << " bonds flipped" << endl;
        }

        if(iter%10==0){
            eners = total_bend_stretch(Pos, mesh, lij_t0, mask_ids, mbrane_para, area_para); 
            fprintf(fp2, "%d %g %g\n", iter, eners.x, eners.y);
            fflush(fp2);
        }

        fprintf(fid , "%d %g", iter, ((float)num_moves/(float)mc_para.one_mc_iter) );
        Ener_t = diag_energies(Et, Pos, Pos_t0, mesh, lij_t0, mbrane_para, area_para, stick_para,
                vol_para,  afm_para,  act_para,  fid );

        outfile_terminal << "iter = " << iter << "; Accepted Moves = " 
            << (double) num_moves*100/mc_para.one_mc_iter << " %;"<<  
            " totalener = "<< mbrane_para.tot_energy[0] << "; volume = " << *mbrane_para.volume <<
            "Area = "<< *mbrane_para.area << endl;

    }

    end_time = MPI_Wtime();
    if(mpi_rank == 0){
        fprintf(stderr, "\n---------\n");
        fprintf(stderr, "Time for %d montecarlo steps = %g min \n", mc_para.tot_mc_iter, (end_time-start_time)/60.0);
        fprintf(stderr, "\n---------\n");
    }

    fclose(fid);
    fclose(fp2);
    free(Pos);
    free(lij_t0);
    free(mesh.node_nbr_list);
    mpi_err = MPI_Finalize();
    outfile_terminal.close();
    return 0;
}


/*     for(iter=0; iter<mbrane_para.N; iter++){ */
/*         printf("%d \n", stick_para.is_attractive[iter]); */
/*     } */

/*     exit(0); */
