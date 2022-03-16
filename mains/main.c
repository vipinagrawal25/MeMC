#include "include/global.h"
#include "include/subroutine.h"
#include <random>
#include <unistd.h>
#include <cmath>
#include <sstream>
#include<iomanip>

template<typename T>
inline string ZeroPadNumber(T num){
    ostringstream ss;
    ss << setw( 5 ) << setfill( '0' ) << (int)num;
    return ss.str();
}

int main(int argc, char *argv[]){

    pid_t pid = getpid();
    cout << "# ID for this process is: " << pid << endl;
    //
    int i, iterations, num_moves;
    double Et[5], Ener_t;
    double vol_sph, e_t, s_t,area_sph;
    POSITION *Pos;
    bool *is_attractive;
    MBRANE_para mbrane;
    MCpara mcpara;
    AFM_para afm;
    MESH mesh;
    POSITION afm_force;
    FILE *fid;
    double *lij_t0
    int *triangles;
    string outfolder, syscmds, log_file, outfile, para_file;
    char log_headers[] = "# iter acceptedmoves total_ener stretch_ener bend_ener stick_ener afm_ener ener_volume  forcex, forcey forcez area nPole_z sPole_z";
    int nPole,sPole;
    int ibydumpskip;
    double Pole_zcoord;
    if(argc!=3){
        printf("\n\n mayday.. requires an argument <parameter file> <output folder>\n\n");
        exit(0);
    }else{
        para_file=argv[1];
        outfolder=argv[2];
    }
    //
    syscmds="mkdir "+outfolder;
    system(syscmds.c_str());
    syscmds="cp "+para_file+" "+outfolder+"/";
    system(syscmds.c_str());
    init_rng();
    // read the input file
    initialize_read_parameters(&mbrane, &afm, &mcpara, para_file);
   /* define all the paras */ 
    mbrane.volume = (double *)calloc(1, sizeof(double)); 
    mbrane.volume[0] = (4./3.)*pi*pow(mbrane.radius,3);
    mbrane.tot_energy = (double *)calloc(1, sizeof(double));
    mbrane.tot_energy[0] = 0e0;
    // allocate arrays
    Pos = (POSITION *)calloc(mbrane.N, sizeof(POSITION));
    mesh.cmlist = (int *)calloc(mbrane.N+1, sizeof(int));
    mesh.node_nbr_list = (int *)calloc(mbrane.num_nbr, sizeof(int)); 
    mesh.bond_nbr_list = (int2 *)calloc(mbrane.num_nbr, sizeof(int2));
    lij_t0 = (double *)calloc(mbrane.num_nbr, sizeof(double));
    triangles = (int *)calloc(mbrane.num_nbr, sizeof(int));
    is_attractive = (bool *)calloc(mbrane.N, sizeof(bool));
    afm.tip_curve = (POSITION *)calloc(afm.N, 3*sizeof(double));
    //
    if(!mcpara.is_restart){
        s_t = afm.sigma; 
        afm.sigma = 0.00;
        e_t = afm.epsilon; 
        afm.epsilon = 0.0;
    }

    if(!mcpara.is_restart){
        hdf5_io_read_config((double *) Pos, (int *) mesh.cmlist,
                (int *) mesh.node_nbr_list, (int2 *) mesh.bond_nbr_list, 
                triangles, "input/input.h5");
        initialize_eval_lij_t0(Pos, mesh, lij_t0, &mbrane);
        identify_attractive_part(Pos, is_attractive, mbrane.N);
        max(&nPole,&Pole_zcoord,Pos,mbrane.N);
        min(&sPole,&Pole_zcoord,Pos,mbrane.N);
    }else{
        hdf5_io_read_config((double *) Pos, (int *) mesh.cmlist,
                (int *) mesh.node_nbr_list, (int2 *) mesh.bond_nbr_list, 
                triangles, "input/input.h5");
        max(&nPole,&Pole_zcoord,Pos,mbrane.N);
        min(&sPole,&Pole_zcoord,Pos,mbrane.N);
        initialize_eval_lij_t0(Pos, mesh, lij_t0, &mbrane);
        identify_attractive_part(Pos, is_attractive, mbrane.N);
        hdf5_io_read_config((double *) Pos, (int *) mes_t.cmlist,
                (int *) mesh.node_nbr_list, (int2 *) mesh.bond_nbr_list,
                triangles, outfolder+"/restart.h5");
    }
    //
    double HH = mbrane.coef_str/(mbrane.av_bond_len*mbrane.av_bond_len);
    double BB = mbrane.coef_bend;
    cout << "# Foppl von Karman (FvK): "
         << 2*HH*mbrane.radius*mbrane.radius/(BB*sqrt(3)) << endl;
    //
    // cout << nPole << endl;
    // cout << Pole_zcoord;
    // exit(1);
    // uncomment these and put N=n*n where no is integer 
    // in afm para to be able to visualize afm tip
    /* initialize_afm_tip(afm); */
    /* sprintf(log_file, "%s/afm_tip.vtk", outfolder); */
    /* visit_vtk_io_afm_tip((double *) afm.tip_curve, */ 
    /*         afm.N, log_file); */
    //
    Et[0] =  stretch_energy_total(Pos, mesh, lij_t0, mbrane);
    Et[1] =  bending_energy_total(Pos, mesh, mbrane);
    Et[2] = lj_bottom_surf_total(Pos, is_attractive, mbrane);
    Et[3] = lj_afm_total(Pos, &afm_force, mbrane, afm);
    //
    volume_area_enclosed_membrane(Pos, triangles, 
            mbrane.num_triangles,&vol_sph,&area_sph);
    double  ini_vol = (4./3.)*pi*pow(mbrane.radius,3);
    Et[4] = mbrane.coef_vol_expansion*(vol_sph/ini_vol - 1e0)*(vol_sph/ini_vol - 1e0);
    Ener_t = Et[0] + Et[1] + Et[2] + Et[3] + Et[4];
    mbrane.tot_energy[0] = Ener_t;
    mbrane.volume[0] = vol_sph;
    //
    log_file=outfolder+"/mc_log";
    fid = fopen(log_file.c_str(), "a");
    if(!mcpara.is_restart)fprintf(fid, "%s\n", log_headers);
    num_moves = 0;
    for(i=0; i < mcpara.tot_mc_iter; i++){
        Et[0] =  stretch_energy_total(Pos, mesh, lij_t0, mbrane);
        Et[1] =  bending_energy_total(Pos, mesh, mbrane);
        Et[2] = lj_bottom_surf_total(Pos, is_attractive, mbrane);
        Et[3] = lj_afm_total(Pos, &afm_force, mbrane, afm);
        volume_area_enclosed_membrane(Pos, triangles, mbrane.num_triangles, &vol_sph, &area_sph);
        Et[4] = mbrane.coef_vol_expansion*(vol_sph/ini_vol - 1e0)*(vol_sph/ini_vol - 1e0);
        // Ener_t = Et[0] + Et[1] + Et[2] + Et[3] + Et[4];
        // mbrane.tot_energy[0] = Ener_t;
        cout << "iter = " << i << "; Accepted Moves = " << (double) num_moves*100/mcpara.one_mc_iter << " %;"<<  
                " totalener = "<< mbrane.tot_energy[0] << "; volume = " << mbrane.volume[0]<< "; area = " << area_sph << endl;
        fprintf(fid, " %d %d %g %g %g %g %g %g %g %g %g %g %g %g\n",
                    i, num_moves, mbrane.tot_energy[0], Et[0], Et[1], Et[2], Et[3], Et[4],
                    afm_force.x, afm_force.y, afm_force.z,area_sph,Pos[nPole].z,Pos[sPole].z);
        fflush(fid);
        if(i%mcpara.dump_skip == 0){
            outfile=outfolder+"/part_"+ ZeroPadNumber(i/mcpara.dump_skip)+".vtk";
            // sprintf(outfile,"%s/part_%05d.vtk",outfolder,(int)i/mcpara.dump_skip);
            visit_vtk_io( (double *) Pos, triangles,
                    mbrane.N, outfile);

            // visit_vtk_io_point_data(is_attractive, mbrane.N,
            //         outfile, "isattr");
            // dump config for restart
            // hdf5_io_dump_restart_config((double *) Pos, (int *) mesh.cmlist,
            //         (int *) mesh.node_nbr_list, (int2 *) mesh.bond_nbr_list,  
            //         triangles, mbrane, outfolder);
            hdf5_io_dump_restart_config((double *) Pos, (int *) mesh.cmlist,
                    (int *) mesh.node_nbr_list, (int2 *) mesh.bond_nbr_list,  
                    triangles, mbrane, "input");
        }
        if(i == 10*mcpara.dump_skip && !mcpara.is_restart ){
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
    free(mesh.bond_nbr_list);
    return 0;
}
