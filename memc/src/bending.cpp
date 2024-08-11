#include "bending.hpp"
#include <fstream>
#define sign(x) ((x > 0) ? 1 : ((x < 0) ? -1 : 0))

extern "C" void BendRead(double *, double *, double *, double *, double *, char *);
int get_nstart(int , int);
/*------------------------*/
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
/*------------------------*/
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
/*------------------------*/
double cotangent(double a, double b, double c){

     /// @brief  a, b, c are the length of the sides of a triangle 
     
     ///  @return   0.25*(a*a+b*b-c*c)/area; where area is the area of triangle
    double s = 0.5*(a+b+c);
    double area = sqrt(s*(s-a)*(s-b)*(s-c));
    double cot_theta=0.25*(a*a+b*b-c*c)/area;
    return cot_theta;
}
/*------------------------*/
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
double BE::bending_energy_ipart(Vec3d *pos, int *node_nbr, int num_nbr, int idx){
    /// @brief Estimate the Bending energy contribution when ith particle 
    /// position changes
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
    double cot_aij[num_nbr],cot_bij[num_nbr],area_ijk[num_nbr];
    double ljksq[num_nbr];
    Vec3d xij[num_nbr],xjk[num_nbr];
    double cot_jdx_k,cot_jdx_kp,cot_kdx,cot_kpdx,area_ijkp;
    Vec3d xik,xikp,xjkp,nhat_local,xijp1;
    int jdx,kdx,kpdx,jdxp1;
    double cot_sum,liksq,likpsq,ljkpsq,lijsq[num_nbr];
    sigma_i = 0e0;
    // store all the lengths
    for (int j = 0; j < num_nbr; ++j){
        jdx = node_nbr[j];
        kdx = node_nbr[(j+1)%num_nbr]; // this is same as kdx
        xij[j] = pos[idx]- pos[jdx];
        xjk[j] = pos[jdx]- pos[kdx];
        lijsq[j] = inner_product(xij[j],xij[j]);
        ljksq[j] = inner_product(xjk[j],xjk[j]);
        area_ijk[j] = 0.5*norm(cross_product(xij[j],xjk[j]));
    }
    // Now compute all the angles
    for (int j = 0; j < num_nbr; ++j){
        liksq=lijsq[(j+1)%num_nbr];
        likpsq=lijsq[(j-1+num_nbr)%num_nbr];
        ljkpsq=ljksq[(j-1+num_nbr)%num_nbr];
        area_ijkp=area_ijk[(j-1+num_nbr)%num_nbr];
        cot_aij[j] = 0.25*(ljkpsq+likpsq-lijsq[j])/area_ijkp;
        cot_bij[j] = 0.25*(ljksq[j]+liksq-lijsq[j])/area_ijk[j];
    }
    for (int j = 0; j < num_nbr; j++){
        cot_jdx_k = cot_aij[(j+1)%num_nbr];
        cot_jdx_kp = cot_bij[(j-1+num_nbr)%num_nbr];
        liksq=lijsq[(j+1)%num_nbr];
        likpsq=lijsq[(j-1+num_nbr)%num_nbr];
        area_ijkp=area_ijk[(j-1+num_nbr)%num_nbr];
        xijp1=xij[(j+1)%num_nbr];
        cot_sum=0.5*(cot_aij[j] + cot_bij[j]);
        cot_times_rij = cot_times_rij + xij[j]*cot_sum;
        sigma_i=sigma_i+voronoi_area(cot_jdx_k,cot_bij[j],liksq,lijsq[j],area_ijk[j]);
        nhat_local=cross_product(xijp1,xij[j]);
        nhat = nhat + nhat_local*(1e0/norm(nhat_local));
    }
    nhat = nhat/norm(nhat);
    lap_bel = cot_times_rij/sigma_i;
    lap_bel_t0 = nhat*curv_t0;
    bend_ener = 0.5*BB*sigma_i*normsq(lap_bel-lap_bel_t0);
    return bend_ener;
}
/*------------------------*/
double BE::bending_energy_ipart_neighbour(Vec3d *pos, MESH_p mesh, int idx){
    /// @brief Estimate the Bending energy contribution from the neighbours 
    /// when ith particle position changes
    /// @param Pos array containing co-ordinates of all the particles
    /// @param mesh mesh related parameters -- connections and neighbours information; 
    /// @param para  Membrane related parameters;
    /// @return Bending Energy contribution from the neighbours of ith particle
    int j;
    int num_nbr_j;
    int nbr, cm_idx_nbr;
    double be;
    be = 0e0;
    for(j = idx*mesh.nghst; j < idx*mesh.nghst + mesh.numnbr[idx]; j++){
        nbr = mesh.node_nbr_list[j];
        num_nbr_j = mesh.numnbr[nbr];
        cm_idx_nbr = nbr*mesh.nghst;
        be += bending_energy_ipart(pos, 
              (int *) mesh.node_nbr_list + cm_idx_nbr,
               num_nbr_j, nbr);
    }
    return be;
}
/*------------------------*/
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
                (int *) (mesh.node_nbr_list + cm_idx), num_nbr, idx);
     }
     return be;
}
/*-------------------------------------------------------------------------------------*/
