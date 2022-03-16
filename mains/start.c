int main(int argc char **argv){
     int i, iterations, num_moves;
     double Ener;
     char *conf;
     POSITION *Pos;
     Neighbours *neib;
     LJpara  para;
     MCpara mcpara;
     FILE *fid;


    /* define all the paras */ 

     para.N  = 5120;
     para.len = 2*pi;
     para.epsilon = 1;
     /* para.sigma = para.len/sqrt((double)para.N); */

     para.sigma = sqrt(8*pi/(2*para.N-4));
     para.r_cut = 4*para.sigma;

     // define the monte carlo parameters
     mcpara.dfac  = 32;
     mcpara.one_mc_iter = 10*para.N;
     mcpara.kBT = 1;
     mcpara.metric = "sph";


     Pos = calloc(para.N, sizeof(POSITION));
     neib = calloc(para.N, sizeof(Neighbours));

     initialize_system(Pos, para, mcpara);


     printf("%lf  %lf \n", para.sigma, para.r_cut);

     make_nlist_pf(Pos, neib,  para, mcpara.metric);


     Ener = pairlj_total_energy_pf(Pos,  para, 
             mcpara.metric);

     iterations = 60000;
     system("touch output/mc_log");
     fid = fopen("output/mc_log", "a");
     for(i=0; i<iterations; i++){
         if(i%1000 == 0){
             fprintf(stderr, " iter, AcceptedMoves, Energy: %d %d %g\n",
                     i, num_moves, Ener);
             dump_config(Pos,  para.len, i, para.N);
         }
         fprintf(fid, " %d %d %g\n",
                 i, num_moves, Ener);
         fflush(fid);
         num_moves = monte_carlo_surf2d(Pos, neib, 
                 para, mcpara);
         make_nlist_pf(Pos, neib,  para, mcpara.metric);

         Ener = pairlj_total_energy(Pos, neib, para, 
                 mcpara.metric);
         /*       /1* "a comment"; *1/ */
     }
     fclose(fid);
     free(Pos);
     free(neib);
     return 0;
 }

