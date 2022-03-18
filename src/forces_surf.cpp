#include "global.h"
#include "subroutine.h"

#define sign(x) ((x > 0) ? 1 : ((x < 0) ? -1 : 0))

double cotangent(Vec3d si, Vec3d sk, Vec3d sj){
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
    double s = 0.5*(a+b+c);
    double area = sqrt(s*(s-a)*(s-b)*(s-c));
    double cot_theta=0.25*(a*a+b*b-c*c)/area;
    return cot_theta;
}

Vec3d determine_xyz_parabola(Vec3d pos, AFM_para afm) {
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
        fprintf(stderr, "X0 small roots may not be correct\n");
    }

    a0 = 1./afm.tip_rad;
    /* a0 = a0*a0; */
    c0 = afm.tip_pos_z;

    a = (2*pow(a0*a0*x0, 2)  + 2*pow(a0*a0*y0,2));
    b = 0e0;
    c = x0*x0 + 2*a0*a0*x0*x0*(c0 - z0);
    d = -x0*x0*x0;

    nroot = cubic_solve(a, b, c, d, roots);

     if(nroot  == 3){
        fprintf(stderr, "Multiple points satisfy distance minimization \n");
    }

    pt_pbola.x = roots[0];
    pt_pbola.y = (y0/x0)*roots[0];
    pt_pbola.z  = (a0*pt_pbola.x)*(a0*pt_pbola.x) + 
        (a0*pt_pbola.y)*(a0*pt_pbola.y) + c0;

    return pt_pbola;

}
double volume_ipart(Vec3d *pos, 
        int *node_nbr, int2* bond_nbr,
        int num_nbr, int idx, MBRANE_para para){

    int i, j, k, kp;
    double volume1, volume2;
    Vec3d area1, area2, rk, ri, rj, rkp;
    Vec3d rij, rijk, rik, rjkp;
    Vec3d rijkp;
    Vec3d rjk, pt; 

    volume1 = 0e0;
    volume2 = 0e0;
    for (i =0; i < num_nbr; i++){
        j = node_nbr[i];
        k  = bond_nbr[i].i1; 
        kp = bond_nbr[i].i2; 
        ri = pos[idx]; rj = pos[j]; 
        rk = pos[k]; rkp = pos[kp];

        rij = Vec3d_add(ri, rj, 1e0);
        rijk = Vec3d_add(rij, rk, 1e0);
        rijkp = Vec3d_add(rij, rkp, 1e0);

        rijk.x = rijk.x/3e0; rijk.y = rijk.y/3e0; rijk.z = rijk.z/3e0;
        rijkp.x = rijkp.x/3e0; 
        rijkp.y = rijkp.y/3e0; rijkp.z = rijkp.z/3e0;

        rij  = Vec3d_add(rj, ri, -1e0);
        rjk  = Vec3d_add(rj , rk, -1e0);
        rjkp = Vec3d_add(rj , rkp, -1e0);


        area1 = cross_product(rij, rjk);
        area2 = cross_product(rij, rjkp);
        double ip1 = inner_product(area1,rijk);
        double ip2 = inner_product(area2,rijkp);

        if(sign(ip1) == 1) volume1 = volume1 + ip1;
        if(sign(ip1) == -1) volume2 = volume2 + ip1;

        if(sign(ip2) == 1) volume1 = volume1 + ip2;
        if(sign(ip2) == -1) volume2 = volume2 + ip2;
    }
    volume1 = volume1/3e0;
    return volume1;
}
double stretch_energy_ipart(Vec3d *pos, 
        int *node_nbr, double *lij_t0,
        int num_nbr, int idx, MBRANE_para para){
    //
    double HH;
    double idx_ener;
    Vec3d rij;
    double mod_rij;
    int i,j;
    //
    idx_ener = 0e0;
    HH = para.coef_str/(para.av_bond_len*para.av_bond_len);
    for (i =0; i < num_nbr; i++){
        j = node_nbr[i];
        rij = Vec3d_add(pos[idx], pos[j], -1e0);
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
double voronoi_area(double cotJ, double cotK, double jsq, double ksq, double area){
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
double bending_energy_ipart(Vec3d *pos, int *node_nbr, int2 *bond_nbr, int num_nbr,
                            int idx, MBRANE_para para, string method){
    double bend_ener,curvature,sigma_i;
    Vec3d cot_times_rij;
    double BB=para.coef_bend;
    double curv_t0 = para.sp_curv;
    // lap_bel:Laplace Beltrami operator for sphere is 2\kappa\nhat
    // nhat is outward normal.
    Vec3d lap_bel,lap_bel_t0,nhat;
    // double curv_t0 = 2e0/para.radius;
    cot_times_rij.x = 0e0;
    cot_times_rij.y = 0e0;
    cot_times_rij.z = 0e0;
    //
    if (method=="voro_negative"){
        double rij_dot_rij;
        int i;
        // double curvature;
        Vec3d cot_theta_rij, rij, rik, rikp;
        //
        sigma_i = 0e0;
        //
        // for (i =0; i < num_nbr; i++){
        //     j = node_nbr[i];
        //     k  = bond_nbr[i].i1; 
        //     kp = bond_nbr[i].i2;
        //     //
        //     cot_sum=0.5*(cotangent(pos[idx],pos[k],pos[j]) +
        //                  cotangent(pos[idx],pos[kp],pos[j]));
        //     rij = Vec3d_add(pos[idx], pos[j], -1e0);
        //     // rik = Vec3d_add(pos[idx], pos[k], -1e0);
        //     // rikp = Vec3d_add(pos[idx], pos[kp], -1e0);
        //     // area = area+sqrt(inner_product(cross_product(rij,rik),cross_product(rij,rik)))+
        //     //             sqrt(inner_product(cross_product(rij,rikp),cross_product(rij,rikp)));
        //     rij_dot_rij = inner_product(rij, rij);
        //     cot_times_rij  = Vec3d_add(cot_times_rij, rij, cot_sum);
        //     // cot_times_rij  = Vec3d_add(cot_times_rij, rij, 1/sqrt(3));
        //     sigma_i = sigma_i + rij_dot_rij*cot_sum;
        // }
        // sigma_i *= 0.25;        
        int jdx,kdx,kpdx;
        Vec3d xij,xik,xjk,xikp,xjkp;
        double cot_kdx,cot_kpdx;
        double cot_sum;
        for (i =0; i < num_nbr; i++){
            jdx = node_nbr[i];
            kdx  = bond_nbr[i].i1; 
            kpdx = bond_nbr[i].i2;
            //
            cot_kdx = cotangent(pos[idx],pos[kdx],pos[jdx]);
            cot_kpdx = cotangent(pos[idx],pos[kpdx],pos[jdx]);
            //
            cot_sum=0.5*(cot_kdx +  cot_kpdx);
            //
            xij = Vec3d_add(pos[idx], pos[jdx], -1e0);
            //
            // rik = Vec3d_add(pos[idx], pos[k], -1e0);
            // rikp = Vec3d_add(pos[idx], pos[kp], -1e0);
            // area = area+sqrt(inner_product(cross_product(rij,rik),cross_product(rij,rik)))+
            //             sqrt(inner_product(cross_product(rij,rikp),cross_product(rij,rikp)));
            rij_dot_rij = inner_product(xij, xij);
            cot_times_rij  = Vec3d_add(cot_times_rij, xij, cot_sum);
            // cot_times_rij  = Vec3d_add(cot_times_rij, rij, 1/sqrt(3));
            sigma_i = sigma_i + 0.25*(cot_sum*rij_dot_rij);
        }
        // is_attractive[idx]=sigma_i<0;
        // if (sigma_i<0){
        // sigma_i=area/3;
        // cout << area/3 << endl;
        // }
    }else if (method=="curvature_unsigned" || method=="curvature"){
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
            kdx  = bond_nbr[j].i1;
            kpdx = bond_nbr[j].i2;
            //
            xij = Vec3d_add(pos[idx], pos[jdx], -1e0);
            xijp1 = Vec3d_add(pos[idx], pos[jdxp1], -1e0);
            xik = Vec3d_add(pos[idx], pos[kdx], -1e0);
            xjk = Vec3d_add(pos[jdx], pos[kdx], -1e0);
            xikp = Vec3d_add(pos[idx], pos[kpdx], -1e0);
            xjkp = Vec3d_add(pos[jdx], pos[kpdx], -1e0);
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
            // cot_idx_kp = 0.25*(lijsq+likpsq-ljkpsq)/area_ijkp;
            // compute sigma_i -> first check whether all angles are acute?
            //
            // cot_jdx_k = cotangent(pos[idx],pos[jdx],pos[kdx]);
            // cot_kdx = cotangent(pos[idx],pos[kdx],pos[jdx]);
            // cot_jdx_kp = cotangent(pos[idx],pos[jdx],pos[kpdx]);
            // cot_kpdx = cotangent(pos[idx],pos[kpdx],pos[jdx]);
            // //
            cot_sum = 0.5*(cot_kdx+cot_kpdx);
            cot_times_rij = Vec3d_add(cot_times_rij, xij, cot_sum);
            // sigma_i = sigma_i + 0.125*( cot_jdx_k*liksq+
            //                             cot_kdx*lijsq+
            //                             cot_jdx_kp*likpsq+
            //                             cot_kpdx*lijsq);
            // // sigma_i=sigma_i+0.125*();
            sigma_i=sigma_i+voronoi_area(cot_jdx_k,cot_kdx,liksq,lijsq,area_ijk);
            sigma_i=sigma_i+voronoi_area(cot_jdx_kp,cot_kpdx,likpsq,lijsq,area_ijkp);
            //
            nhat_local=cross_product(xijp1,xij);
            nhat=Vec3d_add(nhat,nhat_local,1e0/norm(nhat_local));
        }
        nhat = nhat/norm(nhat);
        // nhat=nhat/norm(nhat);
        sigma_i = 0.5*sigma_i;  // as everything is counted twice.
        // double xij[num_nbr],cot_alpha[num_nbr],cot_beta[num_nbr],xij_sq[num_nbr];
        // int jdx,kdx,jpdx,kpdx,jp;
        // double sigma_i=0;
        // // bool obtuse=0;      // if any cotangent is -ve or cot(alpha+beta)<0
        // for (int j = 0; j < num_nbr; ++j){
        //     // jp = (j+1)%(num_nbr);         // things are cyclic
        //     if(cot_alpha[jp] > 0 && cot_beta[j] > 0 && cot_alpha[jp]*cot_beta[j] < 1){
        //         // this means that all the angles are acute.
        //         // add the voronoi area.
        //         sigma_i=sigma_i+0.25*(cot_alpha[jp]*xij_sq[jp] + cot_beta[j]+xij_sq[j]);
        //     }else{ 
        //         // add non-voronoi area
        //         if(cot_alpha[jp]*cot_beta[j] > 1){
        //         }
        //     }
        // }
    }
    if (method=="curvature_unsigned" || method == "voro_negative"){
        curvature = (1e0/sigma_i)*norm(cot_times_rij);
        bend_ener = 0.5*BB*sigma_i*(curvature-curv_t0)*(curvature-curv_t0);
    }else if(method=="curvature"){
        lap_bel = cot_times_rij/sigma_i;
        lap_bel_t0 = nhat*curv_t0;
        bend_ener = 0.5*BB*sigma_i*normsq(lap_bel-lap_bel_t0);
    }
    return bend_ener;
 }

double bending_energy_ipart_neighbour(Vec3d *pos, 
        MESH mesh, int idx, MBRANE_para para){
/*    Evaluates the bending energy of all the neighbours*/ 
/*    of idx particle*/
   int j;
   int num_nbr_j;
   int nbr, cm_idx_nbr;
   double be;

   be = 0e0;

   for (j = mesh.cmlist[idx]; j < mesh.cmlist[idx+1]; j++){
       nbr = mesh.node_nbr_list[j];
       num_nbr_j = mesh.cmlist[nbr+1] -  mesh.cmlist[nbr];
       cm_idx_nbr = mesh.cmlist[nbr];
       be += bending_energy_ipart(pos, 
              (int *) mesh.node_nbr_list + cm_idx_nbr,
              (int2 *) mesh.bond_nbr_list + cm_idx_nbr, num_nbr_j, 
              nbr, para);
   }
   return be;
} 

double bending_energy_total(Vec3d *pos, 
        MESH mesh, MBRANE_para para){
    int idx;
    int num_nbr, cm_idx;
    double be;

    be = 0e0;
    for(idx = 0; idx < para.N; idx++){
        /* idx = 2; */
        num_nbr = mesh.cmlist[idx + 1] - mesh.cmlist[idx];
        cm_idx = mesh.cmlist[idx];


        /* fprintf(stderr, "check neighbour of second: %d %d %d \n", num_nbr, cm_idx); */

        be += bending_energy_ipart(pos, 
                (int *) (mesh.node_nbr_list + cm_idx),
                (int2 *) (mesh.bond_nbr_list + cm_idx), num_nbr, 
                idx, para);
    }
    return be;
}

double stretch_energy_total(Vec3d *pos, 
        MESH mesh, double *lij_t0,
         MBRANE_para para){
    int idx;
    int num_nbr, cm_idx;
    double se;

    se = 0e0;
    for(idx = 0; idx < para.N; idx++){
        /* idx = 2; */
        num_nbr = mesh.cmlist[idx + 1] - mesh.cmlist[idx];
        cm_idx = mesh.cmlist[idx];
        se += stretch_energy_ipart(pos, 
                (int *) (mesh.node_nbr_list + cm_idx),
                (double *) (lij_t0 + cm_idx), num_nbr, 
                idx, para);
        /* printf( "stretch: %lf \n", se); */
    }
    return se*0.5e0;
 }


double lj_rep(double sqdr, double eps){
    double r6;
    r6 = sqdr*sqdr*sqdr;
    return eps*(r6*(r6));
}

double lj_attr(double sqdr, double eps){
    double r6;
    r6 = sqdr*sqdr*sqdr;
    return 4*eps*(r6*(r6-1));
}

double lj_bottom_surface(double zz, 
        bool is_attractive, 
        double sur_pos, double eps, double sigma){

    double inv_sqdz, ds;

    if(is_attractive){
        ds = sur_pos - zz;
        inv_sqdz = (sigma*sigma)/(ds*ds);
        return  lj_attr(inv_sqdz, eps);
    } else {
        return 0e0;
    }
}


double lj_bottom_surf_total(Vec3d *pos, 
        bool *is_attractive, MBRANE_para para){
    int idx;
    double lj_bote;

    lj_bote = 0e0;
    for(idx = 0; idx < para.N; idx++){
        /* idx = 2; */

        lj_bote += lj_bottom_surface(pos[idx].z, is_attractive[idx],
                para.pos_bot_wall, para.epsilon, para.sigma);
    }
    return lj_bote;
}



void identify_attractive_part(Vec3d *pos, 
        bool *is_attractive, int N){

    int i; 
    double theta;
    for(i= 0; i<N; i++){
        theta = pi - acos(pos[i].z);
        is_attractive[i] = theta < pi/6.0;
    }
}

double lj_afm(Vec3d pos, AFM_para afm){
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


double lj_afm_total(Vec3d *pos, 
        Vec3d *afm_force, MBRANE_para para,
        AFM_para afm){
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



void volume_area_enclosed_membrane(Vec3d *pos, 
    int *triangles, int num_triangles,
    double *avolume, double *aarea){
    int i, j, k, it;
    Vec3d area, rk, ri, rj;
    Vec3d rij, rijk, rik;
    Vec3d rjk, pt; 
    double volume,tot_area;

    volume = 0e0;
    tot_area=0e0;
    i = triangles[0];
    for (it = 0; it< 3*num_triangles; it=it+3){
        i = triangles[it];
        j = triangles[it+1];
        k = triangles[it+2];
        ri = pos[i]; rj = pos[j]; rk = pos[k];
        /* ri = Vec3d_add(pos[i], pt, -1); */
        /* rj = Vec3d_add(pos[j], pt, -1); */
        /* rk = Vec3d_add(pos[k], pt, -1); */
        
        rij = Vec3d_add(ri, rj, 1e0);
        rijk = Vec3d_add(rij, rk, 1e0);

        rijk.x = rijk.x/3e0; rijk.y = rijk.y/3e0; rijk.z = rijk.z/3e0;

        rij = Vec3d_add(rj, ri, -1e0);
        rik = Vec3d_add(rk , ri, -1e0);

        area = cross_product(rij, rik);
        volume = volume + 0.5*fabs(area.x*rijk.x + area.y*rijk.y + area.z*rijk.z);
        tot_area += 0.5*norm(area);
    }
    *aarea = tot_area;
    *avolume = volume/3e0;
};

void identify_obtuse(Vec3d *pos, int *triangles,
       double *obtuse,  int N){
    double piby2;
    int i, j, k, it;
    double a_ij_ik, a_ji_jk, a_ki_kj;
    Vec3d ri, rj, rk;
    Vec3d rij, rik, rkj;
    double inp_r1, inp_r2, inp_r3;

    piby2 = 0.5*pi;

    for (it = 0; it< 3*N; it=it+3){
        i = triangles[it];
        j = triangles[it+1];
        k = triangles[it+2];

        obtuse[it/3] = 0e0;
        ri = pos[i]; rj = pos[j]; rk = pos[k];

        rij = Vec3d_add(ri, rj, -1e0);
        inp_r1 = inner_product(rij, rij);

        rik = Vec3d_add(ri , rk, -1e0);
        inp_r2 = inner_product(rik, rik);

        rkj = Vec3d_add(rk , rj, -1e0);
        inp_r3 = inner_product(rkj, rkj);

        a_ij_ik = (inner_product(rij,rik))/(sqrt(inp_r1*inp_r2));
        a_ij_ik = acos(a_ij_ik);

        a_ji_jk = (inner_product(rij,rkj))/(sqrt(inp_r1*inp_r3));
        a_ji_jk = acos(a_ji_jk);

        a_ki_kj = -(inner_product(rik,rkj))/(sqrt(inp_r2*inp_r3));
        a_ki_kj = acos(a_ki_kj);

       bool logic = ((a_ij_ik > piby2) ||
                (a_ji_jk > piby2) || 
                (a_ki_kj > piby2));

       /* double sum_ang = a_ij_ik + a_ji_jk + a_ki_kj; */

       /* printf( "%lf %lf %lf %lf %lf\n", a_ij_ik, */
               /* a_ji_jk, a_ki_kj, sum_ang, pi); */ 
      
        if(logic) obtuse[it/3] = 1e0;
    }
}
