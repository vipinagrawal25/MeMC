#include "../includes/global.h"
#include "../includes/subroutine.h"



int main(int argc, char *argv[]){
    //
    int i, iterations, num_moves;
    double Et[5], Ener_t;
    double vol_sph, e_t, s_t;
    Vec3d *Pos;
    bool *is_attractive;
    MBRANE_para mbrane;
    MCpara mcpara;
    AFM_para afm;
    MESH mesh;
    Vec3d afm_force;
    FILE *fid;
    double *lij_t0;
    char *outfolder, *syscmds, *log_file, *outfile, *para_file;

    outfolder = (char *)malloc(128*sizeof(char));
    syscmds = (char *)malloc(128*sizeof(char));
    log_file = (char *)malloc(128*sizeof(char));
    outfile = (char *)malloc(128*sizeof(char));
    para_file = (char *)malloc(128*sizeof(char));
    char log_headers[] = "# iter acceptedmoves total_ener stretch_ener bend_ener stick_ener afm_ener ener_volume";

    if(argc!=3){
        printf("\n\n mayday.. requires an argument <parameter file> <output folder>\n\n");
        exit(0);
    }else{
        para_file=argv[1];
        outfolder=argv[2];
    }
    //
    sprintf(syscmds,"mkdir %s",outfolder);
    system(syscmds);
    sprintf(syscmds,"%s %s %s%s", (char *)"cp", para_file,outfolder,(char *)"/");
    system(syscmds);
    init_rng(23397);
    // read the input file
    init_read_parameters(&mbrane, &afm, &mcpara, para_file);
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
    s_t = afm.sigma; 
    afm.sigma = 0.00;
    e_t = afm.epsilon; 
    afm.epsilon = 0.0;
    hdf5_io_read_pos( (double *)Pos,  (char *) "conf/dmemc_pos.h5");
    hdf5_io_read_mesh((int *) mesh.cmlist,
            (int *) mesh.node_nbr_list,  (char *) "conf/dmemc_conf.h5");
    init_eval_lij_t0(Pos, mesh, lij_t0, &mbrane);
    identify_attractive_part(Pos, is_attractive, mbrane.theta, mbrane.N);
    //
    
    Et[0] = stretch_energy_total(Pos, mesh, lij_t0, mbrane);
    Et[1] = bending_energy_total(Pos, mesh, mbrane);
    Et[2] = lj_bottom_surf_total(Pos, is_attractive, mbrane);
    Et[3] = lj_afm_total(Pos, &afm_force, mbrane, afm);
    //
    vol_sph = volume_total(Pos, mesh, mbrane);

    double  ini_vol = (4./3.)*pi*pow(mbrane.radius,3);
    Et[4] = mbrane.coef_vol_expansion*(vol_sph/ini_vol - 1e0)*(vol_sph/ini_vol - 1e0);
    Ener_t = Et[0] + Et[1] + Et[2] + Et[3] + Et[4];
    mbrane.tot_energy[0] = Ener_t;
    mbrane.volume[0] = vol_sph;
    //
    sprintf(log_file,"%s/mc_log",outfolder);
    fid = fopen(log_file, "a");
    fprintf(fid, "%s\n", log_headers);
    num_moves = 0;
    for(i=0; i < mcpara.tot_mc_iter; i++){
        Et[0] =  stretch_energy_total(Pos, mesh, lij_t0, mbrane);
        Et[1] =  bending_energy_total(Pos, mesh, mbrane);
        Et[2] = lj_bottom_surf_total(Pos, is_attractive, mbrane);
        Et[3] = lj_afm_total(Pos, &afm_force, mbrane, afm);
        vol_sph = volume_total(Pos, mesh, mbrane);
        Et[4] = mbrane.coef_vol_expansion*(vol_sph/ini_vol - 1e0)*(vol_sph/ini_vol - 1e0);
        fprintf(stderr, "iter :: %d AcceptedMoves% :: %4.2f total energy :: %g volume :: %g \n", i, 100*(double)num_moves/mcpara.one_mc_iter, mbrane.tot_energy[0], vol_sph);

        fprintf(fid, " %d %d %g %g %g %g %g %g\n",
                    i, num_moves, mbrane.tot_energy[0], Et[0], Et[1], Et[2], Et[3], Et[4]);
        fflush(fid);
        if(i%mcpara.dump_skip == 0){

            sprintf(outfile,"%s/snap_%04d.h5",outfolder,(int)(i/mcpara.dump_skip));
            hdf5_io_write_pos((double*) Pos, 3*mbrane.N, outfile);
        }
        if(i == 10*mcpara.dump_skip){
            afm.sigma = s_t;
            afm.epsilon = e_t;
            e_t = lj_afm_total(Pos, &afm_force, mbrane, afm);
            mbrane.tot_energy[0] += e_t;
        }
        num_moves = monte_carlo_3d(Pos, mesh, lij_t0, is_attractive, 
                mbrane, mcpara, afm);
    }
    fclose(fid);
    free(Pos);
    free(lij_t0);
    free(mesh.node_nbr_list);
    return 0;
}
