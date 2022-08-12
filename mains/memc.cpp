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
    MESH mesh, mes_t;
    Vec3d afm_force,spring_force[2];
    FILE *fid;
    double *lij_t0;
    double Pole_zcoord;
    string outfolder,syscmds, para_file, log_file, outfile;
    // outfolder = (char *)malloc(128*sizeof(char));
    // syscmds = (char *)malloc(128*sizeof(char));
    // log_file = (char *)malloc(128*sizeof(char));
    // outfile = (char *)malloc(128*sizeof(char));
    // para_file = (char *)malloc(128*sizeof(char));
    double YY=mbrane.YY;
    double BB = mbrane.coef_bend;
    char log_headers[] = "# iter acceptedmoves total_ener stretch_ener bend_ener stick_ener afm_ener ener_volume";
    SPRING_para spring;    
    if(argc!=3){
        printf("\n\n mayday.. requires an argument <parameter file> <output folder>\n\n");
        exit(0);
    }else{
        para_file=argv[1];
        outfolder=argv[2];
    }
    /**** create folder and copy parameter file *****/
    syscmds="mkdir outfolder";
    if(system(syscmds.c_str()) != 0) fprintf(stderr, "failure in creating folder\n");
    syscmds="cp "+para_file+" "+outfolder+"/";
    if(system(syscmds.c_str()) != 0) fprintf(stderr, "failure in copying parafile\n");
    write_param(outfolder + "/para.out",mbrane,mcpara,spring);
    /*************************************************/
    init_rng(23397);
    // read the input file
    init_read_parameters(&mbrane, &afm, &mcpara, &spring, para_file);
   /* define all the paras */ 
    mbrane.volume = (double *)calloc(1, sizeof(double)); 
    mbrane.volume[0] = (4./3.)*pi*pow(mbrane.radius,3);
    mbrane.tot_energy = (double *)calloc(1, sizeof(double));
    mbrane.tot_energy[0] = 0e0;
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
        hdf5_io_read_pos( (double *)Pos,  "conf/dmemc_pos.h5");
        hdf5_io_read_mesh((int *) mesh.cmlist,
                (int *) mesh.node_nbr_list, "conf/dmemc_conf.h5");
        init_eval_lij_t0(Pos, mesh, lij_t0, &mbrane, &spring);
        identify_attractive_part(Pos, is_attractive, mbrane.theta, mbrane.N);
        max(&mesh.nPole,&Pole_zcoord,Pos,mbrane.N);
        min(&mesh.sPole,&Pole_zcoord,Pos,mbrane.N);
    }else{
        hdf5_io_read_pos( (double *)Pos,   outfolder+"/dmemc_pos.h5");
        hdf5_io_read_mesh((int *) mesh.cmlist,
                (int *) mesh.node_nbr_list,  outfolder+"/dmemc_conf.h5");
        max(&mesh.nPole,&Pole_zcoord,Pos,mbrane.N);
        min(&mesh.sPole,&Pole_zcoord,Pos,mbrane.N);
        init_eval_lij_t0(Pos, mesh, lij_t0, &mbrane, &spring);
        identify_attractive_part(Pos, is_attractive, mbrane.theta, mbrane.N);
        hdf5_io_read_pos( (double *)Pos,  outfolder+"/restart.h5");
    }
    /***** Copy initial things to the folder *****/
    if (!mcpara.is_restart){
        syscmds="cp conf/dmemc_pos.h5 "+outfolder+"/";
        system(syscmds.c_str());
        syscmds="cp conf/dmemc_conf.h5 "+outfolder+"/";
        system(syscmds.c_str());
    }
    /*****  initialize energy values *****/
    Et[0] = stretch_energy_total(Pos, mesh, lij_t0, mbrane);
    Et[1] = bending_energy_total(Pos, mesh, mbrane);
    exit(1);
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
    /************************************/
    cout << "# Foppl von Karman (FvK): " 
         << YY*mbrane.radius*mbrane.radius/BB << endl;
    //
    // sprintf(log_file,"%s/mc_log",outfolder);
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
        cout << "iter = " << i << "; Accepted Moves = " << (double) num_moves*100/mcpara.one_mc_iter << " %;"<<  
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
                mbrane, mcpara, afm, spring);
    }
    fclose(fid);
    free(Pos);
    free(lij_t0);
    free(mesh.node_nbr_list);
    return 0;
}