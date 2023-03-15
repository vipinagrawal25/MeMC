#include "global.h"
#include "subroutine.h"
#include <random>
#include <unistd.h>

int main(int argc, char **argv){
     pid_t pid = getpid();
     cout << "# ID for this process is: " << pid << endl;
     //
     int i,  num_moves;
     double Ener;
     Vec2d *Pos;
     Nbh_list *neib;
     LJ_p  para;
     MC_p  mcpara;
     FILE *fid;
     string outfolder,outfile,syscmds;
     char *metric;
     char log_file[128];
     //
     metric = (char *) malloc(128*sizeof(char));
     // outfolder = (char *) malloc(128*sizeof(char));
     // outfile = (char *) malloc(128*sizeof(char));
     if(argc!=5){
         fprintf(stderr, "\n\n Requires argument <Number of points> <output folder>\n\n");
         exit(0);
     }else{
         para.N=atoi(argv[1]);
         metric=argv[2];
         outfolder=argv[3];
         mcpara.tot_mc_iter = atoi(argv[4]);
         fprintf(stderr, "Monte Carlo of %d particles on %s grid..\n",
                 para.N, metric);
         
     }
     //
     syscmds="mkdir "+outfolder;
    // sprintf(syscmds, "%s %s","mkdir ",outfolder);
    if(system(syscmds.c_str()) != 0) fprintf(stderr, "failure in creating folder");
    init_rng(23177);
    /* define all the paras */ 
     para.len = 2*pi;
     para.epsilon = 1;
     para.bdry_condt = 3;
     // 0 for channel; 1 for frame; default is periodic;
     if(strcmp(metric,"cart")==0){
         para.sigma = para.len/sqrt((double)para.N);
     }else if(strcmp(metric,"sph")==0){
         para.sigma = sqrt(8*pi/(2*para.N-4));
     }

     para.r_cut = 4*para.sigma;
     mcpara.dfac  = 16;
     mcpara.one_mc_iter = 2*para.N;
     mcpara.kBT = 1;
     mcpara.dump_skip = 100;

     Pos = (Vec2d *)calloc(para.N, sizeof(Vec2d));
     neib = (Nbh_list *) calloc(para.N, sizeof(Nbh_list));

     init_system_random_pos(Pos, para.len, para.N, metric, para.bdry_condt );

     fprintf(stderr, "For N = %d, sigma = %lf, and epsilon=%lf \n", para.N,
             para.sigma, para.epsilon);

     make_nlist(Pos, neib,  para, metric);


    Ener = pairlj_total_energy(Pos, neib,  para, metric);

    sprintf(log_file, "%s%s",outfolder.c_str(),"/mc_log");
    fid = fopen(log_file, "a");

    for(i=0; i<mcpara.tot_mc_iter; i++){
        if(i%mcpara.dump_skip == 0){
          outfile=outfolder+"/snap_"+ZeroPadNumber(i/mcpara.dump_skip)+".h5";
            // sprintf(outfile, "%s%s%04d%s",outfolder,"/snap_",(int)(i/mcpara.dump_skip),".h5");
            hdf5_io_write_pos((double *) Pos, 2*para.N, outfile);
        }
        fprintf(fid, " %d %d %g\n", i, num_moves, Ener);
        fflush(fid);
        num_moves = monte_carlo_surf2d(Pos, neib, 
                para, mcpara, metric);
        Ener = pairlj_total_energy(Pos, neib, para, 
                metric);
        fprintf(stderr, " iter = %d percentage of AcceptedMoves = %4.2f Energy = %g\n",
                i, 100*(double)num_moves/mcpara.one_mc_iter, Ener);
        make_nlist(Pos, neib,  para, metric);
    }
    fclose(fid);
    free(Pos);
    free(neib);
    return 0;
}

