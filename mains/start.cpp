#include "global.h"
#include "subroutine.h"
#include <random>
#include <unistd.h>



int main(int argc, char **argv){
     int i,  num_moves;
     double Ener;
     Vec2d *Pos;
     Nbh_list *neib;
     LJpara  para;
     MCpara mcpara;
     FILE *fid;
     char *metric, *outfolder, *outfile;
     char syscmds[128], log_file[128];

     metric = (char *) malloc(128*sizeof(char));
     outfolder = (char *) malloc(128*sizeof(char));
     outfile = (char *) malloc(128*sizeof(char));
     if(argc!=4){
         fprintf(stderr, "\n\n Requiresi argument <Number of points> <output folder>\n\n");
         exit(0);
     }else{
         para.N=atoi(argv[1]);
         metric=argv[2];
         outfolder=argv[3];
         fprintf(stderr, "Monte Carlo of %d particles on %s grid..\n",
                 para.N, metric);
         
     }
     //
    sprintf(syscmds, "%s %s","mkdir ",outfolder);
    system(syscmds);
   

    init_rng(23077);
    /* define all the paras */ 

     para.len = 2*pi;
     para.epsilon = 1;
     if(strcmp(metric,"cart")==0){
         para.sigma = para.len/sqrt((double)para.N);
     }else if(strcmp(metric,"sph")==0){
         para.sigma = sqrt(8*pi/(2*para.N-4));
     }

     para.r_cut = 4*para.sigma;
     mcpara.dfac  = 32;
     mcpara.one_mc_iter = 2*para.N;
     mcpara.kBT = 1;
     mcpara.tot_mc_iter = 60000;
     mcpara.dump_skip = 100;

     Pos = (Vec2d *)calloc(para.N, sizeof(Vec2d));
     neib = (Nbh_list *) calloc(para.N, sizeof(Nbh_list));

     init_system_random_pos(Pos, para.len, para.N, metric);

     fprintf(stderr, "For N = %d, sigma = %lf, and epsilon=%lf \n", para.N,
             para.sigma, para.epsilon);

     make_nlist(Pos, neib,  para, metric);


    Ener = pairlj_total_energy(Pos, neib,  para, metric);

    sprintf(log_file, "%s%s",outfolder,"/mc_log");
    fid = fopen(log_file, "a");

    for(i=0; i<mcpara.tot_mc_iter; i++){
        if(i%mcpara.dump_skip == 0){
            sprintf(outfile, "%s%s%04d%s",outfolder,"/part_pos",(int)(i/mcpara.dump_skip),".bin");
            io_dump_config((double *) Pos, 2*para.N, outfile);
        }
        fprintf(fid, " %d %d %g\n", i, num_moves, Ener);
        fflush(fid);
        num_moves = monte_carlo_surf2d(Pos, neib, 
                para, mcpara, metric);
        Ener = pairlj_total_energy(Pos, neib, para, 
                metric);
        fprintf(stderr, " iter :: %d AcceptedMoves :: %d Energy ::  %g\n",
                i, num_moves, Ener);
        make_nlist(Pos, neib,  para, metric);
    }
    fclose(fid);
    free(Pos);
    free(neib);
    return 0;
}

