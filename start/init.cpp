#include "global.h"
#include "subroutine.h"
std::mt19937 rng2;

void init_system_random_pos(Vec2d *Pos,  double len_x, double len_y, 
        int N, char *metric, int bdry_condt ){
    /// @brief Initializes the points on surface of sphere or flat plane
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param metric Topology of the surface "cart" for flat plane "sph" for
    /// sphere
    ///  @param len length of the domain;  
    ///  @param N number of points;
    bool is_sph, is_cart;
    double dl_x, dl_y;

    // this should be calculated or passed parameters
    int n_ghost, n_ghost_x, n_ghost_y;
    // remove it once debugged 

    n_ghost = (int)sqrt(N);
    //distribute particles according to side length
    n_ghost_x = (int)2*n_ghost*len_x/(len_x+len_y);
    n_ghost_y = (int)2*n_ghost*len_y/(len_x+len_y);
    //n_ghost_x=n_ghost;
    //n_ghost_y=n_ghost;
    dl_x = (len_x/(n_ghost_x));
    dl_y = (len_y/(n_ghost_y)); //original: len_y/n_ghost. accordingly changed the scale

    std::uniform_real_distribution<> rand_x(dl_x, len_x-dl_x);
    std::uniform_real_distribution<> rand_y(dl_y, len_y-dl_y);

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
            n_ghost = 2*n_ghost_x;
            for(int i=0; i<n_ghost/2; i++){
                    Pos[i].x = i*dl_x;
                    Pos[i].y = 0.0; 
            }
            for(int i=n_ghost/2; i<n_ghost; i++){
                Pos[i].x = (i - n_ghost/2 + 0.5)*dl_x;
                /* Pos[i].x = (i - n_ghost/2)*(2*(len + 0.5)/n_ghost);  bp*/
                Pos[i].y = len_y; 
            }
            for(int i=n_ghost; i<N; i++){
                Pos[i].x = rand_x(rng2);
                Pos[i].y = rand_y(rng2);
            }
            break;
      case 1:
            // This is a frame ;
            // n_ghost has to be even for logic to work;
            n_ghost = 2*n_ghost_x + 2*n_ghost_y;
            for(int i=0; i<n_ghost_x; i++){
                Pos[i].x = (i+0.5)*dl_x;
                Pos[i].y = 0.0; 
            }
            for(int i=n_ghost_x; i<2*n_ghost_x; i++){
                Pos[i].x = (i - n_ghost_x + 0.5)*dl_x;
                Pos[i].y = len_y; 
            }
            for(int i=2*n_ghost_x; i<2*n_ghost_x+n_ghost_y; i++){
                Pos[i].x = 0.0;
                Pos[i].y = (i - 2*n_ghost_x + 0.5)*dl_y; 
            }

            for(int i=2*n_ghost_x+n_ghost_y; i<n_ghost; i++){
                Pos[i].x = len_x;
                Pos[i].y = (i +0.5 -2*n_ghost_x-n_ghost_y)*dl_y; 
            }
            for(int i=n_ghost; i<N; i++){
                Pos[i].x = rand_x(rng2);
                Pos[i].y = rand_y(rng2);
            }

            break;
       
            default:
                for(int i=0; i<N; i++){
                    Pos[i].x = rand_x(rng2);
                    Pos[i].y = rand_y(rng2);
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


