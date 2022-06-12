#include "global.h"
#include "subroutine.h"


double cal_length(double x1 , double x2, 
        double y1, double y2, double len, char
        *metric){

	 ///  @brief Calculate the length between points x1,y1 and x2,y2 
  	 ///
	 ///  @param metric Topology of the surface "cart" for flat plane "sph" for
     /// sphere
	 ///  @param x1  Coordinate x1 if metric is cart; Theta_1 if metric is sphere;
	 ///  @param x2  Coordinate x2 if metric is cart; Theta_2 if metric is sphere;
	 ///  @param y1  Coordinate y1 if metric is cart; Phi_1 if metric is sphere;
	 ///  @param y2  Coordinate y2 if metric is cart; Phi_2 if metric is sphere;
	 ///  @param len length of the domain;
     
	 ///  @return   (x2-x1)^2 + (y2-y1)^2 if metric is cart;
     ///   (sin(x2)cos(y2) - sin(x1)*cos(y1))^2 +  
     ///   (sin(x2)sin(y2) - sin(x1)*sin(y1))^2 +  
     ///   (cos(x2) - cos(x1))^2  if metric is sphere
	 /// 

    bool is_sph, is_cart;
    double dx, dy, dz, ds;
    double spx1, spx2, spy1, spy2, spz1, spz2;

    is_sph = false;
    is_cart = false;
    if(strcmp(metric, "sph") == 0){
        is_sph = true;
    }
    if(strcmp(metric, "cart") == 0){
        is_cart = true;
    }

    if (is_cart){
        dx = x2 - x1;
        dy = y2 - y1;
        if(dx >  len*0.5) dx = dx - len; 
        if(dy >  len*0.5) dy = dy - len;
        if(dx < -len*0.5) dx = dx + len;
        if(dy < -len*0.5) dy = dy + len;
        ds = dx*dx + dy*dy;
    }
    if (is_sph){
        spx1 = sin(x1)*cos(y1);
        spy1 = sin(x1)*sin(y1);
        spz1 = cos(x1);

        spx2 = sin(x2)*cos(y2);
        spy2 = sin(x2)*sin(y2);
        spz2 = cos(x2);

        dx = spx2 - spx1;
        dy = spy2 - spy1;
        dz = spz2 - spz1;

        ds = dx*dx + dy*dy + dz*dz;

    }

    return ds;

}


void make_nlist(Vec2d *Pos, Nbh_list *neib,
        LJpara para, char *metric){

	 ///  @brief Makes a list of neighbouring particles for every particle  
  	 ///
	 ///  @param Pos array containing co-ordinates of all the particles
	 ///  @param neib array containing array of neigbours for each particle  
	 ///  @param para  Parameters of lennard jones potential
	 ///  @param metric Topology of the surface "cart" for flat plane "sph" for
     /// sphere
     ///  @details calculates distance between the particle and stores the index of
     /// particle lying within a circle of certain radius in neib

    double  new_rc;


    for(int i=0;i < para.N; i++)
        neib[i].cnt = 0;


    new_rc = para.r_cut*para.r_cut;

    for(int i = 0; i < para.N; i++){
        for(int j = i+1; j < para.N; j++){
            int m = i;
            int n = j;
            if(len_check(Pos[m], Pos[n],  para.len, new_rc, metric)){
                neib[m].list[neib[m].cnt] = n;
                neib[n].list[neib[n].cnt] = m;
                neib[m].cnt += 1;
                neib[n].cnt += 1;
            }
        }

    }
}


bool len_check(Vec2d s1, Vec2d s2,
        double len, double new_rc, char *metric){

	 ///  @brief checks the distance between  s1 and s2 
     ///  @param s1 coordinate of first particle
     ///  @param s2 coordinate of second particle
     /// @param new_rc square of the cutoff length
	 ///  @param metric Topology of the surface "cart" for flat plane "sph" for
     /// sphere
     /// @param len length of the domain
  	 /// @return  
     /// True if (s2 - s1)^2 is smaller than ds, False otherwise

    double ds;
    bool var = false;

    ds = cal_length(s1.x , s2.x, 
            s1.y, s2.y, len, metric);
    var = (ds < new_rc);
    return var;
}

double pairlj_ipart_energy(Vec2d *Pos, int *n_list,
        int ni, int i_p, LJpara para, char *metric){
     /// @brief Estimate the energy of ith particle
	 ///  @param Pos array containing co-ordinates of all the particles
	 ///  @param i_p index of ith particle;
	 ///  @param n_list array containing index of particles neibouring i_p; 
	 ///  @param ni number of neigbours;
	 ///  @param para  Parameters of lennard jones potential
	 ///  @param metric Topology of the surface "cart" for flat plane "sph" for
     /// sphere
     /// @return Energy of the ith particle.
     ///  @details see https://en.wikipedia.org/wiki/Lennard-Jones_potential


    double Sq_dr2;
    int j, k;
    double Sq_dr1,r2,r6;
    double loc_ener;
    double r2_cut;
    double inv_sig_ma, eps, ds;


    loc_ener = 0e0;
    for(j=0; j < ni; j++){
        k = n_list[j];
        if( k != i_p ){
            ds = cal_length(Pos[i_p].x , Pos[k].x, 
                    Pos[i_p].y , Pos[k].y, para.len, metric);
            Sq_dr1 = ds;
            inv_sig_ma = 1e0/(para.sigma*para.sigma);
            eps = para.epsilon;
            r2_cut = para.r_cut*para.r_cut;
            Sq_dr2 = Sq_dr1*inv_sig_ma;
            if(Sq_dr1 < r2_cut){	
                r2 = 1.0/Sq_dr2;
                r6  = r2*r2*r2; 
                loc_ener += eps*r6*(r6);
            }
        }
    }
    return loc_ener;
}


double pairlj_total_energy(Vec2d *Pos, Nbh_list *neib,
        LJpara para, char *metric){
   /// @brief Estimate the total energy of the system
	 ///  @param Pos array containing co-ordinates of all the particles
	 ///  @param neib array containing list of neibouring particles to each particle; 
	 ///  @param para  Parameters of lennard jones potential
	 ///  @param metric Topology of the surface "cart" for flat plane "sph" for
     /// sphere
     /// @return Total energy of the system.
     ///  @details see https://en.wikipedia.org/wiki/Lennard-Jones_potential


    int i,ni;
    double j_pe, pe;

    pe = 0e0;
    for(i=0;i<para.N;i++){
        ni = neib[i].cnt;
        j_pe = pairlj_ipart_energy(Pos, neib[i].list,
                ni, i, para, metric);
        pe += j_pe;
    }
    return 0.5*pe;
}
