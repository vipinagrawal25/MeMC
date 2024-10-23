#include "bending.hpp"
#include "multicomp.hpp"
#include <fstream>

#define sign(x) ((x > 0) ? 1 : ((x < 0) ? -1 : 0))

extern "C" void BendRead(double *, double *, double *, double *, char *);
int get_nstart(int , int);
/*------------------------*/
BE::BE(const MESH_p& mesh, std::string fname){
    char tmp_fname[128];
    string parafile, outfile;

    parafile = fname+"/para_file.in";
    sprintf(tmp_fname, "%s", parafile.c_str());
    BendRead(&bend1, &bend2, &spC1, &spC2, tmp_fname);
    init_coefbend(mesh.compA, mesh.N);
    spcurv=spC1;

    if (mesh.ncomp==1){
        bend2=bend1;
        spC2=spC1;
    }

    ofstream out_;
    out_.open( fname+"/bendpara.out");
    out_<< "# =========== bending parameters ==========" << endl
      << " N " << mesh.N << endl
      << " bend1 = " << bend1 << endl
      << " bend2 = " << bend2 << endl
      << " spC1 " << spC1 << endl
      << " spC2 " << spC2 << endl;
    out_.close();
}
/*-------------------------------------------------*/
void BE::init_coefbend(int *lipA, int N){
    for (int i = 0; i < N; ++i) {
        // std::cout << lipA[i] << std::endl;
        if (lipA[i]) coef_bend.push_back(bend2);
        else coef_bend.push_back(bend1);
    }
}
/*-------------------------------------------------*/
double BE::bending_energy_ipart(Vec3d *pos, int *node_nbr, int num_nbr, int idx,
            int bdry_type, double lenth, int edge){
    /// @brief Estimate the Bending energy contribution when ith particle 
    /// position changes
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param idx index of ith particle;
    ///  @param node_nbr nearest neigbours of idx; 
    ///  @param num_nbr number of nearest neigbours of idx; 
    ///  @param para  Membrane related parameters;
    /// @todo try openMP Pragmas;
    /// @return Bending Energy contribution when ith particle is displaced.
    double bend_ener;
    Vec3d cot_times_rij;
    Vec3d lap_bel,lap_bel_t0,nhat;
    double cot_aij[num_nbr],cot_bij[num_nbr],area_ijk[num_nbr];
    double ljksq[num_nbr], lijsq[num_nbr];
    Vec3d xij[num_nbr],xjk;
    double cot_jdx_k,cot_jdx_km,cot_kdx,cot_kmdx,area_ijkm;
    Vec3d xik,xikm,xjkm,nhat_local,xijp1;
    int jdx,kdx,kmdx,jdxp1;
    double cot_sum,liksq,likmsq,ljkmsq;
    double sigma_i = 0e0;
    // store all the lengths
    if (bdry_type == 1 || idx>edge){
        for (int j = 0; j < num_nbr; ++j){
            jdx = node_nbr[j];
            kdx = node_nbr[(j+1)%num_nbr]; // this is same as kdx
            xij[j]=pos[idx]-pos[jdx];
            xjk = pos[jdx]-pos[kdx];
            lijsq[j] = inner_product(xij[j],xij[j]);
            ljksq[j] = inner_product(xjk,xjk);
            area_ijk[j] = 0.5*norm(cross_product(xij[j],xjk));
        }
    }else{
        for (int j = 0; j < num_nbr; ++j){
            jdx = node_nbr[j];
            kdx = node_nbr[(j+1)%num_nbr]; // this is same as kdx
            xij[j] = diff_pbc(pos[idx],pos[jdx],lenth);
            xjk =  diff_pbc(pos[jdx],pos[kdx], lenth);
            lijsq[j] = inner_product(xij[j],xij[j]);
            ljksq[j] = inner_product(xjk,xjk);
            area_ijk[j] = 0.5*norm(cross_product(xij[j],xjk));
        }
    }
    // Now compute all the angles
    for (int j = 0; j < num_nbr; ++j){
        liksq=lijsq[(j+1)%num_nbr];
        likmsq=lijsq[(j-1+num_nbr)%num_nbr];
        ljkmsq=ljksq[(j-1+num_nbr)%num_nbr];
        area_ijkm=area_ijk[(j-1+num_nbr)%num_nbr];
        cot_aij[j] = 0.25*(ljkmsq+likmsq-lijsq[j])/area_ijkm;
        cot_bij[j] = 0.25*(ljksq[j]+liksq-lijsq[j])/area_ijk[j];
    }
    for (int j = 0; j < num_nbr; j++){
        cot_jdx_k = cot_aij[(j+1)%num_nbr];
        cot_jdx_km = cot_bij[(j-1+num_nbr)%num_nbr];
        liksq=lijsq[(j+1)%num_nbr];
        likmsq=lijsq[(j-1+num_nbr)%num_nbr];
        area_ijkm=area_ijk[(j-1+num_nbr)%num_nbr];
        xijp1=xij[(j+1)%num_nbr];
        cot_sum=0.5*(cot_aij[j] + cot_bij[j]);
        cot_times_rij = cot_times_rij + xij[j]*cot_sum;
        sigma_i=sigma_i+voronoi_area(cot_jdx_k,cot_bij[j],liksq,lijsq[j],area_ijk[j]);
        nhat_local=cross_product(xijp1,xij[j]);
        nhat = nhat + nhat_local*(1e0/norm(nhat_local));
    }
    nhat = nhat/norm(nhat);
    lap_bel = cot_times_rij/sigma_i;
    lap_bel_t0 = nhat*spcurv;
    bend_ener = 0.5*coef_bend[idx]*sigma_i*normsq(lap_bel-lap_bel_t0);
    return bend_ener;
}
/*--------------------------------------------------------------------------*/
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
            num_nbr_j, nbr, mesh.bdry_type, mesh.boxlen, mesh.edge);
    }
    return be;
}
/*------------------------*/
double BE::bending_energy_total(Vec3d *pos, MESH_p mesh){
    /// @brief Estimate the total Bending energy
    ///  @param Pos array containing co-ordinates of all the particles
    /// @param mesh mesh related parameters -- connections and 
    /// neighbours information;
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
        be += bending_energy_ipart(pos, (int *) (mesh.node_nbr_list + cm_idx),
                num_nbr, idx, mesh.bdry_type, mesh.boxlen, mesh.edge);
    }
    return be;
}
/*------------------------------------------------------------------------------*/
void BE::exchange(int idx1, int idx2){
   double temp = coef_bend[idx1];
   coef_bend[idx1] = coef_bend[idx2];
   coef_bend[idx2] = temp;
}
// /*-------------------------------------------------------------------------------------*/
// double BE::Gaussian_energy_ipart(Vec3d *pos, int *node_nbr, int num_nbr, int idx,
//     int bdry_type, double lenth, int edge){

// }