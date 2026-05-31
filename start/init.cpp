#include "global.h"
#include "subroutine.h"
std::mt19937 rng2;

void init_system_random_pos(Vec2d *Pos,  double len, 
        int N, char *metric, int bdry_condt ){
    /// @brief Initializes the points on surface of sphere or flat plane
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param metric Topology of the surface "cart" for flat plane "sph" for
    /// sphere
    ///  @param len length of the domain;
    ///  @param N number of points;
    bool is_sph, is_cart;
    double dl;

    // this should be calculated or passed parameters
    int n_ghost;
    // remove it once debugged 

    n_ghost = (int) sqrt(N);
    dl = (len/n_ghost);

    std::uniform_real_distribution<> rand_real(dl, len-dl);

    is_sph = false;
    is_cart = false;
    if(strcmp(metric, "sph") == 0){
        is_sph = true;
    }
    if(strcmp(metric, "cart") == 0){
        is_cart = true;
    }
    if(!is_cart && !is_sph){
        fprintf(stderr, "Unknown metric, input should be: \n");
        fprintf(stderr, "a) cart for a 2d cartesian monte carlo \n");
        fprintf(stderr, "b) sph  for a monte carlo on the surface of a sphere\n");
        exit(0);
    }

     
    if(is_cart){
        switch (bdry_condt) {
            case 0:
            // This is a channel ; read bdry_condt as 0 
            // n_ghost has to be even for logic to work;
            n_ghost = 2*n_ghost;
            for(int i=0; i<n_ghost/2; i++){
                    Pos[i].x = i*dl;
                    Pos[i].y = 0.0; 
            }
            for(int i=n_ghost/2; i<n_ghost; i++){
                Pos[i].x = (i - n_ghost/2 + 0.5)*dl;
                /* Pos[i].x = (i - n_ghost/2)*(2*(len + 0.5)/n_ghost);  bp*/
                Pos[i].y = len; 
            }
            for(int i=n_ghost; i<N; i++){
                Pos[i].x = rand_real(rng2);
                Pos[i].y = rand_real(rng2);
            }
            break;
      case 1:
            // This is a frame ;
            // n_ghost has to be even for logic to work;
            n_ghost = 4*n_ghost;
            for(int i=0; i<n_ghost/4; i++){
                Pos[i].x = (i+0.5)*dl;
                Pos[i].y = 0.0; 
            }
            for(int i=n_ghost/4; i<n_ghost/2; i++){
                Pos[i].x = (i - n_ghost/4 + 0.5)*dl;
                Pos[i].y = len; 
            }
            for(int i=n_ghost/2; i<3*n_ghost/4; i++){
                Pos[i].x = 0.0;
                Pos[i].y = (i - n_ghost/2 + 0.5)*dl; 
            }

            for(int i=3*n_ghost/4; i<n_ghost; i++){
                Pos[i].x = len;
                Pos[i].y = (i +0.5 - 3*n_ghost/4)*dl; 
            }
            for(int i=n_ghost; i<N; i++){
                Pos[i].x = rand_real(rng2);
                Pos[i].y = rand_real(rng2);
            }

            break;
       
            default:
                for(int i=0; i<N; i++){
                    Pos[i].x = rand_real(rng2);
                    Pos[i].y = rand_real(rng2);
                }
        }
    }
    if(is_sph){
        Pos[0].x = 0;
        Pos[0].y = 0;
        Pos[1].x = pi;
        Pos[1].y = 0;
        for(int i=2; i<N; i++){
            Pos[i].x = acos(2*drand48() - 1); 
            Pos[i].y = 2*pi*drand48();
        }
        Pos[2].x = acos(2*drand48() - 1); 
        Pos[2].y = 2*pi*drand48();

    }
    /* for(int i=0; i<N; i++){ */
    /*     printf("%lf %lf \n", Pos[i].x); */
    /*     Pos[i].x = drand48()*len; */
    /*     Pos[i].y = drand48()*len; */
    /* } */
}
/*--------------------------------------------------------------------------------*/


