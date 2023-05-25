#include "Vector.h"
#include "global.h"
#include "subroutine.h"
#define sign(x) ((x > 0) ? 1 : ((x < 0) ? -1 : 0))
double cotangent(Vec3d si, Vec3d sk, Vec3d sj){
     ///
     ///  @param si  coordinate of ith point
     ///  @param sk  coordinate of kth point
     ///  @param sj  coordinate of jth point 
     
     ///  @return   ({si-sk}.{sj-sk})/sqrt(({si-sk}x{sj-sk})^2)
     /// angle between vector si and sj
     ///
    Vec3d drik, drjk, cross;
    double cot_theta;  
    double inner_prod;
    //
    drik = Vec3d_add(si , sk, -1e0);
    drjk = Vec3d_add(sj , sk, -1e0);
    cross = cross_product(drik, drjk);
    inner_prod = inner_product(drik, drjk);
    cot_theta = inner_prod/sqrt(inner_product(cross,cross));
    //
    return cot_theta;
}
//
double cotangent(double a, double b, double c){

     /// @brief  a, b, c are the length of the sides of a triangle 
     
     ///  @return   0.25*(a*a+b*b-c*c)/area; where area is the area of triangle
    double s = 0.5*(a+b+c);
    double area = sqrt(s*(s-a)*(s-b)*(s-c));
    double cot_theta=0.25*(a*a+b*b-c*c)/area;
    return cot_theta;
}

Vec3d determine_xyz_parabola(Vec3d pos, AFM_p afm) {
    ///
     ///  @param pos  position of a point in membrane
     ///  @param afm  parameters of AFM tip
     //
     /// 
     ///  @return  position of point on afm tip nearest to pos 
     ///


    int nroot;
    double a, b, c, d;
    double a0, c0;
    double roots[6];
    Vec3d pt_pbola;
    double x0, y0, z0;
    // a coef of x^3
    // b coef of x^2
    // c coef of x^1
    // d coef of x^0

    x0 = pos.x;
    y0 = pos.y;
    z0 = pos.z;

    // some exception messages
    if(fabs(x0) < 1e-15){
        if(fabs(y0)< 1e-15 & fabs(z0)< 1e-15 ){
            fprintf(stderr, "x0 = %g, y0 = %g,  z0 = %g; Points representing exosome very small...\n", x0, y0, z0);
            fprintf(stderr, "See section Checking execution status in https://github.com/vipinagrawal25/MeMC/blob/main/README.md \n");
        }
        x0 = x0 + 2e-12*(2*(drand48() - 1));
    }

    a0 = 1./afm.tip_rad;
    /* a0 = a0*a0; */
    c0 = afm.tip_pos_z;

    a = (2*pow(a0*a0*x0, 2)  + 2*pow(a0*a0*y0,2));
    b = 0e0;
    c = x0*x0 + 2*a0*a0*x0*x0*(c0 - z0);
    d = -x0*x0*x0;

    nroot = cubic_solve(a, b, c, d, roots);

/*      if(nroot  == 3){ */
/*         fprintf(stderr, "Multiple points satisfy distance minimization \n"); */
/*     } */

    pt_pbola.x = roots[0];
    pt_pbola.y = (y0/x0)*roots[0];
    pt_pbola.z  = (a0*pt_pbola.x)*(a0*pt_pbola.x) + 
        (a0*pt_pbola.y)*(a0*pt_pbola.y) + c0;

    return pt_pbola;

}

double area_ipart(Vec3d *pos, int *node_nbr,
        int num_nbr, int idx){
     /// @brief Estimate the area substended by the ith particle
     ///  @param Pos array containing co-ordinates of all the particles
     ///  @param idx index of ith particle;
     ///  @param node_nbr nearest neigbours of idx; 
     ///  @param num_nbr number of nearest neigbours of idx; 
     ///  @param para  Membrane related parameters;
     /// @return Volume substended by ith particle.

    int i, j, k;
    Vec3d  rk, ri, rj, rkp;
    Vec3d rij, ip1, rik, rjkp;
    Vec3d rjk, pt; 
    double area1, area2;

    area1 = 0e0;
    ri = pos[idx];
    for (i =0; i < num_nbr; i++){
        j = node_nbr[i];
        k=node_nbr[(i+1)%num_nbr];
        ri = pos[idx]; rj = pos[j]; 
        rk = pos[k];
        rij  = ri - rj;
        rjk  = ri - rk;
        ip1 = cross_product(rjk, rij);
        area1 = area1 + 0.5*norm(ip1);
    }
    return area1;
}


double volume_ipart(Vec3d *pos, int *node_nbr,
        int num_nbr, int idx){
     /// @brief Estimate the volume substended by voronoi area of the ith particle
     ///  @param Pos array containing co-ordinates of all the particles
     ///  @param idx index of ith particle;
     ///  @param node_nbr nearest neigbours of idx; 
     ///  @param num_nbr number of nearest neigbours of idx; 
     ///  @param para  Membrane related parameters;
     /// @return Volume substended by ith particle.


    int i, j, k;
    double volume1;
    Vec3d area1, area2, rk, ri, rj, rkp;
    Vec3d rij, rijk, rik, rjkp;
    Vec3d rijkp;
    Vec3d rjk, pt; 

    volume1 = 0e0;
    ri = pos[idx];
    for (i =0; i < num_nbr; i++){
        j = node_nbr[i];
        k=node_nbr[(i+1)%num_nbr];
        ri = pos[idx]; rj = pos[j]; 
        rk = pos[k];
        rij = Vec3d_add(ri, rj, 1e0);
        rijk = Vec3d_add(rij, rk, 1e0);
        rijk.x = rijk.x/3e0; rijk.y = rijk.y/3e0; rijk.z = rijk.z/3e0;
        rij  = Vec3d_add(ri, rj, -1e0);
        rjk  = Vec3d_add(ri , rk, -1e0);
        area1 = cross_product(rjk, rij);
        double ip1 = 0.5*inner_product(area1,rijk);
        volume1 = volume1 + abs(ip1);
    }
    volume1 = volume1/3e0;
    return volume1;
}

 double stretch_energy_ipart(Vec3d *pos,
         int *node_nbr, double *lij_t0,
         int num_nbr, int idx, AREA_p para){

    /// @brief Estimate the Stretching energy contribution when ith particle position changes
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param idx index of ith particle;
    ///  @param node_nbr nearest neigbours of idx; 
    /// @param lij_t0 initial distance between points of membrane
    ///  @param num_nbr number of nearest neigbours of idx; 
    ///  @param para  Membrane related parameters;
    /// @return Stretching Energy contribution when ith particle is displaced 


    //
    double HH;
    double idx_ener;
    Vec3d rij;
    double mod_rij, avg_lij;
    int i,j;
    //
    idx_ener = 0e0;
    // HH = para.coef_str/(para.av_bond_len*para.av_bond_len);
    HH = para.YY*sqrt(3)/2;
    for (i =0; i < num_nbr; i++){
        j = node_nbr[i];
        rij = pos[idx] - pos[j];
        mod_rij = sqrt(inner_product(rij, rij));
        //
        idx_ener = idx_ener + (mod_rij - lij_t0[i])*(mod_rij - lij_t0[i]);
    }
    //
    return 0.5*idx_ener*HH;
}

//
//Q. Given two cotangent angles, it returns either the area due to perpendicular bisector,
// or the barycenter.
double voronoi_area(double cotJ, double cotK, 
        double jsq, double ksq, double area){
    /// @brief Estimate the area of the voronoi cell. If I J K are the nodes of
    /// triangle
    ///  @param cotJ angle at node j (see paper/paper.pdf)
    ///  @param cotK angle at node k (see paper/paper.pdf)
    ///  @param jsqr square of the length of bond i-k
    /// @param ksq square of the length of bond i-j  
    ///  @param area area of the triangle 
    /// @return  Given two cotangent angles, it returns either the area due to perpendicular bisector,
/// or the barycenter.


    double sigma;
    if (cotJ>0 && cotK>0){
        if (cotJ*cotK<1){
            // all angles are acute;
            sigma = 0.125*(cotJ*jsq+cotK*ksq);
        }else{
            sigma = 0.5*area;
        }
    }else{
       sigma = 0.25*area;
   }
    return sigma;
}
//
//Would like this routine to return the area and Energy
Vec2d bending_energy_ipart(Vec3d *pos, int *node_nbr, int num_nbr,
                            int idx, MBRANE_p para){

    /// @brief Estimate the Bending energy contribution when ith particle position changes
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param idx index of ith particle;
    ///  @param node_nbr nearest neigbours of idx; 
    ///  @param num_nbr number of nearest neigbours of idx; 
    ///  @param para  Membrane related parameters;
    /// @todo try openMP Pragmas;
    /// @return Bending Energy contribution when ith particle is displaced 

    double bend_ener,sigma_i;
    Vec3d cot_times_rij;
    double BB=para.coef_bend;
    double curv_t0 = para.sp_curv;
    // lap_bel:Laplace Beltrami operator for sphere is 2\kappa\nhat
    // nhat is outward normal.
    Vec3d lap_bel,lap_bel_t0,nhat;
    Vec2d bend_ener_area;
    double vor_area;
    //
    double cot_jdx_k,cot_jdx_kp,cot_kdx,cot_kpdx;
    double area_ijk,area_ijkp;
    double lijsq,liksq,ljksq,likpsq,ljkpsq;
    // Vec3d cot_times_rij;
    Vec3d xij,xik,xjk,xikp,xjkp,nhat_local,xijp1;
    int jdx,kdx,kpdx,jdxp1;
    double cot_sum;
    sigma_i = 0e0;
    for (int j = 0; j < num_nbr; j++){
        jdx = node_nbr[j];
        jdxp1=node_nbr[(j+1)%num_nbr];
        //
        /* kdx  = bond_nbr[j].i1; */
        kdx = node_nbr[(j+1+num_nbr)%num_nbr];
        /* kpdx = bond_nbr[j].i2; */
        kpdx = node_nbr[(j-1+num_nbr)%num_nbr];
        //
        xij = pos[idx]- pos[jdx];
        xijp1 = pos[idx]- pos[jdxp1];
        xik = pos[idx]- pos[kdx];
        xjk = pos[jdx]- pos[kdx]; 
        xikp = pos[idx]- pos[kpdx];
        xjkp = pos[jdx]- pos[kpdx];
        //
       //
        lijsq = inner_product(xij,xij);
        liksq = inner_product(xik,xik);
        ljksq = inner_product(xjk,xjk);
        likpsq = inner_product(xikp,xikp);
        ljkpsq = inner_product(xjkp,xjkp);
        //
        area_ijk = 0.5*norm(cross_product(xij,xjk));
        area_ijkp = 0.5*norm(cross_product(xij,xjkp));
        //
        cot_jdx_k = 0.25*(lijsq+ljksq-liksq)/area_ijk;
        cot_kdx = 0.25*(ljksq+liksq-lijsq)/area_ijk;
        // cot_idx_k = 0.25*(lijsq+liksq-ljksq)/area_ijk;
        //
        cot_jdx_kp = 0.25*(lijsq+ljkpsq-likpsq)/area_ijkp;
        cot_kpdx =  0.25*(ljkpsq+likpsq-lijsq)/area_ijkp;
        cot_sum = 0.5*(cot_kdx+cot_kpdx);
        cot_times_rij = Vec3d_add(cot_times_rij, xij, cot_sum);
        sigma_i=sigma_i+voronoi_area(cot_jdx_k,cot_kdx,liksq,lijsq,area_ijk);
        sigma_i=sigma_i+voronoi_area(cot_jdx_kp,cot_kpdx,likpsq,lijsq,area_ijkp);
        //
        nhat_local = cross_product(xijp1,xij);
        nhat = Vec3d_add(nhat,nhat_local,1e0/norm(nhat_local));
    }
    nhat = nhat/norm(nhat);
    sigma_i = 0.5*sigma_i;  // as everything is counted twice.
    lap_bel = cot_times_rij/sigma_i;
    lap_bel_t0 = nhat*curv_t0;
    bend_ener = 0.5*BB*sigma_i*normsq(lap_bel-lap_bel_t0);
    bend_ener_area.x = bend_ener;
    bend_ener_area.y = sigma_i;
    return bend_ener_area;
}

double bending_energy_ipart_neighbour(Vec3d *pos, 
        MESH_p mesh, int idx, MBRANE_p para){

    /// @brief Estimate the Bending energy contribution from the neighbours when ith particle position changes
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param mesh mesh related parameters -- connections and neighbours
    /// information; 
    ///  @param para  Membrane related parameters;
    /// @return Bending Energy contribution from the neighbours of ith particle 

   int j;
   int num_nbr_j;
   int nbr, cm_idx_nbr;
   double be;
   Vec2d  tmp_be_ar;

   be = 0e0;

   for (j = idx*mesh.nghst; j < idx*mesh.nghst + mesh.numnbr[idx]; j++){
       nbr = mesh.node_nbr_list[j];
       num_nbr_j = mesh.numnbr[nbr];
       cm_idx_nbr = nbr*mesh.nghst;
       tmp_be_ar = bending_energy_ipart(pos, 
              (int *) mesh.node_nbr_list + cm_idx_nbr,
               num_nbr_j, nbr, para);
       be += tmp_be_ar.x;

   
   }
   return be;
} 


 double area_total(Vec3d *pos, MESH_p mesh,
         MBRANE_p para){
     ///  @param Pos array containing co-ordinates of all the particles
     ///  @param mesh mesh related parameters -- connections and neighbours
     /// information;
     ///  @param para  Membrane related parameters;
     /// @return Total volume of the shell


     int idx, st_idx;
     int num_nbr, cm_idx;
     double area;

     st_idx = get_nstart(para.N, para.bdry_type);
     area = 0e0;
     for(idx = 0; idx < para.N; idx++){
         /* idx = 2; */
         cm_idx = idx*mesh.nghst;
         num_nbr = mesh.numnbr[idx];

       area += area_ipart(pos,
                 (int *) (mesh.node_nbr_list + cm_idx),
                  num_nbr, idx);
     }
     return area/3e0;
}


 double volume_total(Vec3d *pos, MESH_p mesh,
         MBRANE_p para){
     /// @brief Estimate the total volume of the shell
     ///  @param Pos array containing co-ordinates of all the particles
     ///  @param mesh mesh related parameters -- connections and neighbours
     /// information;
     ///  @param para  Membrane related parameters;
     /// @return Total volume of the shell


     int idx;
     int num_nbr, cm_idx;
     double vol;

     vol = 0e0;
     for(idx = 0; idx < para.N; idx++){
         /* idx = 2; */
         cm_idx = idx*mesh.nghst;
         num_nbr = mesh.numnbr[idx];

       vol += volume_ipart(pos,
                 (int *) (mesh.node_nbr_list + cm_idx),
                  num_nbr, idx);
     }
     return vol/3e0;
}


 double bending_energy_total(Vec3d *pos, MESH_p mesh, 
         MBRANE_p para){

     /// @brief Estimate the total Bending energy
     ///  @param Pos array containing co-ordinates of all the particles
     ///  @param mesh mesh related parameters -- connections and neighbours
     /// information;
     ///  @param para  Membrane related parameters;
     /// @return Total Bending energy


     int idx, st_idx;
     int num_nbr, cm_idx;
     double be;
     Vec2d be_ar;

     be = 0e0;
     st_idx = get_nstart(para.N, para.bdry_type);
     for(idx = st_idx; idx < para.N; idx++){
         /* idx = 2; */

         cm_idx = idx*mesh.nghst;
         num_nbr = mesh.numnbr[idx];

         be_ar = bending_energy_ipart(pos,
                 (int *) (mesh.node_nbr_list + cm_idx),
                  num_nbr, idx, para);
         be += be_ar.x;
     }
     return be;
}


 double stretch_energy_total(Vec3d *pos,
       MESH_p mesh, double *lij_t0, MBRANE_p para, AREA_p area_p){

    /// @brief Estimate the total Stretching energy  
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param mesh mesh related parameters -- connections and neighbours
    /// information; 
    /// @param lij_t0 initial distance between points of membrane
    ///  @param para  Membrane related parameters;
    /// @return Total Stretching energy 


    int idx, st_idx;
    int num_nbr, cm_idx;
    double se;

    st_idx = get_nstart(para.N, para.bdry_type);
    se = 0e0;
    for(idx = st_idx; idx < para.N; idx++){
        /* idx = 2; */
        num_nbr = mesh.numnbr[idx];
        cm_idx = idx*mesh.nghst;

        se += stretch_energy_ipart(pos,
                 (int *) (mesh.node_nbr_list + cm_idx),
                 (double *) (lij_t0 + cm_idx), num_nbr,
                 idx, area_p);

        /* printf( "stretch: %lf \n", se); */
    }
    return se*0.5e0;
 }


double lj_rep(double sqdr, double eps){

    /// @param sqdr square of the distance between two points
    /// @param eps coefficient of the potential 
    /// @return repulsive part of the Lennard-Jones potential.
    ///  @details see https://en.wikipedia.org/wiki/Lennard-Jones_potential


    double r6;
    r6 = sqdr*sqdr*sqdr;
    return eps*(r6*(r6));
}

double lj_attr(double sqdr, double eps){

    /// @param sqdr square of the distance between two points
    /// @param eps coefficient of the potential 
    /// @return Energy evaluated using Lennard-Jones potential.
    ///  @details see https://en.wikipedia.org/wiki/Lennard-Jones_potential


    double r6;
    r6 = sqdr*sqdr*sqdr;
    return 4*eps*(r6*(r6-1));
}

double stick_bottom_surface(Vec3d pt, Vec3d pt_t0, STICK_p st_p){

    double inv_sqdz, ds;

        ds = norm(pt_t0 - pt);
        inv_sqdz = (st_p.sigma*st_p.sigma)/(ds*ds);
        if(st_p.is_pot_harmonic){
            return  st_p.sigma*ds*ds;
        }else{
            return  lj_attr(inv_sqdz, st_p.epsilon);
        }
}



double lj_bottom_surface(double zz, 
        bool is_attractive, 
        double sur_pos, double eps, double sigma){

    /// @brief Sticking of the bottom surface using LJ potential 
    /// @param zz z-coordinate of a point in shell
    /// @param eps coefficient of the potential 
    /// @param sigma sigma in the expression of LJ potential 
    /// @param sur_pos position of the bottom wall 
    /// @param is_attractive true if the zz sees the  bottom wall 
    /// @return Energy evaluated using Lennard-Jones potential.
    /// @details see https://en.wikipedia.org/wiki/Lennard-Jones_potential


    double inv_sqdz, ds;

    if(is_attractive){
        ds = sur_pos - zz;
        inv_sqdz = (sigma*sigma)/(ds*ds);
        return  lj_attr(inv_sqdz, eps);
    } else {
        return 0e0;
    }
}


double lj_bottom_surf_total(Vec3d *pos, Vec3d *pos_t0,
         MBRANE_p para, STICK_p st_p){

    /// @brief Estimate the total Sticking  
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param is_attractive true for all the particles which sees bottom wall 
    ///  @param para  Membrane related parameters;
    /// @return Total Energy contribution from Sticking to bottom surface

    int idx;
    double lj_bote;

    lj_bote = 0e0;
    for(idx = 0; idx < para.N; idx++){
        /* idx = 2; */

        if(st_p.is_attractive[idx])
            lj_bote += stick_bottom_surface(pos[idx], pos_t0[idx], st_p);

            /* lj_bote += lj_bottom_surface(pos[idx].z, st_p.is_attractive[idx], */
                    /* st_p.pos_bot_wall, st_p.epsilon, st_p.sigma); */
    }
    return lj_bote;
}


double lj_afm(Vec3d pos, AFM_p afm){

    /// @brief Energy contribution from AFM tip to ith point in membrane 
    ///  @param Pos co-ordinates of the ith particle 
    ///  @param afm afm related parameters
    /// @return Energy contribution when form tip to the ith particle
    //
    double ener_afm, ds;
    double ds_sig_inv;
    Vec3d dr, pt_pbola;
    ener_afm = 0e0;
    
    pt_pbola = determine_xyz_parabola(pos, afm);
    if(fabs(afm.tip_pos_z - pt_pbola.z) < 4*afm.sigma) {
        dr = Vec3d_add(pt_pbola, pos, -1); 
        ds = (inner_product(dr,dr));
        ds_sig_inv = (afm.sigma*afm.sigma)/ds;
        ener_afm = lj_rep(ds_sig_inv, afm.epsilon);
    }
   return ener_afm;

};

double lj_afm_total(Vec3d *pos, Vec3d *afm_force,
        MBRANE_p para, AFM_p afm){

    /// @brief Energy contribution from AFM tip to all the point in membrane 
    ///  @param Pos co-ordinates of the ith particle 
    ///  @param para Membrane related parameters
    ///  @param afm_force Total force exerted by afm tip on shell.
    ///  @param afm afm related parameters
    /// @return  Total Energy contribution form the tip to the particle
    //
    int idx;
    double  ds;
    double lj_afm_e;
    double lj_afm_t;
    Vec3d f_t, pt_pbola, dr;

    lj_afm_e = 0e0;
    f_t.x = 0; f_t.y = 0; f_t.z = 0;
    for(idx = 0; idx < para.N; idx++){

        lj_afm_t = lj_afm(pos[idx], afm);
        pt_pbola = determine_xyz_parabola(pos[idx], afm);
        dr = Vec3d_add(pt_pbola, pos[idx], -1); 
        ds = (inner_product(dr,dr));
        f_t.x += 12*lj_afm_t*dr.x/ds;
        f_t.y += 12*lj_afm_t*dr.y/ds;
        f_t.z += 12*lj_afm_t*dr.z/ds;

        lj_afm_e  += lj_afm_t; 

        /* lj_afm_e  = determine_xyz_parabola(pos[idx], afm); */
    }
    *afm_force  = f_t;
    return lj_afm_e;
}

double spring_tot_frame_ener(Vec3d *pos, MBRANE_p mbrane,
        SHEAR_p shear){

    double kk = shear.constant;
    double shift_y = shear.slope;
    double ener_spr, xmin;
    int st_idx = get_nstart(mbrane.N, mbrane.bdry_type);
    int i;
    ener_spr = 0e0;
    for(i=0; i<st_idx; i++){
        double xmin = pos[i].x + shift_y*pos[i].y;
        ener_spr += -kk*pow((pos[i].x - xmin),2)/2;
    }
    return ener_spr;
}
/*
double spring_tot_energy_force(Vec3d *Pos, Vec3d *spring_force, 
                               MESH_p mesh, SPRING_p spring){
    double kk = spring.constant;
    double nZeq = spring.nPole_eq_z;
    double sZeq = spring.sPole_eq_z;
    double ener_spr = kk*pow((Pos[mesh.nPole].z-nZeq),2)/2 +
                      kk*pow((Pos[mesh.sPole].z-sZeq),2)/2 ;
    spring_force[0].z = kk*(nZeq-Pos[mesh.nPole].z);
    spring_force[1].z = kk*(sZeq-Pos[mesh.sPole].z);
    return ener_spr;
}
*/
//
//
double vol_energy_change(MBRANE_p mbrane, VOL_p volp, double dvol){
    double KAPPA = volp.coef_vol_expansion;
    double de_vol=0.0;
    double ini_vol = (4./3.)*pi*pow(mbrane.radius,3);
    if (fabs(KAPPA)>1e-16){
        de_vol = (2*dvol/(ini_vol*ini_vol))*(mbrane.volume[0]  - ini_vol)
            + (dvol/ini_vol)*(dvol/ini_vol);
       de_vol = KAPPA*de_vol;
    }
    return de_vol;
}
//
double PV_change(double pressure, double dvol){ 
    return pressure*dvol;
}
/*
double spring_energy(Vec3d pos, int idx, MESH_p mesh, SPRING_p spring){
    if (spring.icompute==0) return 0;
    double ener_spr=0e0;
    double kk=spring.constant;
    double nZeq = spring.nPole_eq_z;
    double sZeq = spring.sPole_eq_z;
    if (mesh.nPole==idx){
        ener_spr=kk*pow((pos.z-nZeq),2)/2;
    }
    if (mesh.sPole==idx){
        ener_spr=kk*pow((pos.z-sZeq),2)/2;
    }
    return ener_spr;
}
*/
//

double frame_spring_energy(Vec3d pos, Vec3d pos_t0, SHEAR_p shear){
    double ener_spr=0e0;
    double kk=shear.constant;
    double shift_y = shear.slope;

    double xmin = pos_t0.x + shift_y*(pos_t0.y - 3.14159); 
    ener_spr = kk*pow((pos.x - xmin),2)/2;

    return ener_spr;
}

