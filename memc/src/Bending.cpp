#include "Bending.hpp"
#include <fstream>
#define sign(x) ((x > 0) ? 1 : ((x < 0) ? -1 : 0))

extern "C" void BendRead(double *, double *, double *, double *, double *, char *);
int get_nstart(int , int);

int BE::initBE(int N, std::string fname){
  char tmp_fname[128];
  string parafile, outfile;

  parafile = fname+"/para_file.in";
  sprintf(tmp_fname, "%s", parafile.c_str());
  BendRead(&coef_bend, &minC, &maxC, &theta, &spcurv, tmp_fname);

  ofstream out_;
  out_.open( fname+"/bendpara.out");
  out_<< "# =========== bending parameters ==========" << endl
      << " N " << N << endl
      << " coef_bend = " << coef_bend << endl
      << " minC " << minC << endl
      << " maxC " << maxC << endl
      << " theta " << theta << endl
      << " spcurv " << spcurv << endl;
  out_.close();
   return 0;
}


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
    drik = si - sk; 
    drjk = sj - sk; 
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
/*-------------------------------------------------*/
double BE::bending_energy_ipart(Vec3d *pos, int *node_nbr, int num_nbr,
                            int idx){
    /// @brief Estimate the Bending energy contribution when ith particle position changes
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param idx index of ith particle;
    ///  @param node_nbr nearest neigbours of idx; 
    ///  @param num_nbr number of nearest neigbours of idx; 
    ///  @param para  Membrane related parameters;
    /// @todo try openMP Pragmas;
    /// @return Bending Energy contribution when ith particle is displaced.
    double bend_ener,sigma_i;
    Vec3d cot_times_rij;
    double BB=coef_bend;
    double curv_t0 = spcurv;
    Vec3d lap_bel,lap_bel_t0,nhat;
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
        cot_times_rij = cot_times_rij + xij*cot_sum;
        sigma_i=sigma_i+voronoi_area(cot_jdx_k,cot_kdx,liksq,lijsq,area_ijk);
        sigma_i=sigma_i+voronoi_area(cot_jdx_kp,cot_kpdx,likpsq,lijsq,area_ijkp);
        //
        nhat_local=cross_product(xijp1,xij);
        nhat=nhat + nhat_local*(1./norm(nhat_local));
    }
    nhat = nhat/norm(nhat);
    sigma_i = 0.5*sigma_i;  // as everything is counted twice.
    lap_bel = cot_times_rij/sigma_i;
    lap_bel_t0 = nhat*curv_t0;
    bend_ener = 0.5*BB*sigma_i*normsq(lap_bel-lap_bel_t0);
    // cout << BB << " " << bend_ener << " " << normsq(lap_bel)<< endl;
    return bend_ener;
}
//
double BE::bending_energy_ipart_neighbour(Vec3d *pos, 
        MESH_p mesh, int idx){
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
   be = 0e0;
   //
   for (j = idx*mesh.nghst; j < idx*mesh.nghst + mesh.numnbr[idx]; j++){
       nbr = mesh.node_nbr_list[j];
       num_nbr_j = mesh.numnbr[nbr];
       cm_idx_nbr = nbr*mesh.nghst;
       be += bending_energy_ipart(pos, 
              (int *) mesh.node_nbr_list + cm_idx_nbr,
               num_nbr_j, nbr);
   }
   return be;
}

double BE::bending_energy_total(Vec3d *pos, MESH_p mesh){
     /// @brief Estimate the total Bending energy
     ///  @param Pos array containing co-ordinates of all the particles
     ///  @param mesh mesh related parameters -- connections and neighbours
     /// information;
     ///  @param para  Membrane related parameters;
     /// @return Total Bending energy
     int idx, st_idx;
     int num_nbr, cm_idx;
     double be;
     be = 0e0;
     //
     st_idx = get_nstart(mesh.N, mesh.bdry_type);

     for(idx = st_idx; idx < mesh.N; idx++){
         /* idx = 2; */
         cm_idx = idx*mesh.nghst;
         num_nbr = mesh.numnbr[idx];

         be += bending_energy_ipart(pos,
                 (int *) (mesh.node_nbr_list + cm_idx),
                  num_nbr, idx);
     }
     return be;
}
/*-------------------------------------------------------------------------------------*/
