#include "bending.hpp"
#include "multicomp.hpp"
#include <fstream>
#include <cmath>
#include <iterator>

#define sign(x) ((x > 0) ? 1 : ((x < 0) ? -1 : 0))

extern "C" void BendRead(double *, double *, double *, double *, bool*, char *);
int get_nstart(int , int);
/*------------------------*/
BE::BE(const MESH_p& mesh, std::string fname){
    char tmp_fname[128];
    string parafile, outfile;

    parafile = fname+"/para_file.in";
    sprintf(tmp_fname, "%s", parafile.c_str());
    BendRead(&bend1, &bend2, &spC1, &spC2, &iGauss, tmp_fname);
    
    spcurv=spC1;
    
    if (mesh.ncomp==1){
        bend2=bend1;
        spC2=spC1;
        iGauss=false;
    }

    if (bend1 != bend2 && mesh.ncomp>1) multicomp=true;
    ghost=mesh.nghst;   // not recommended.

    ofstream out_;
    out_.open( fname+"/bendpara.out");
    out_<< "# =========== bending parameters ==========" << endl
      << " N " << mesh.N << endl
      << " bend1 = " << bend1 << endl
      << " bend2 = " << bend2 << endl
      << " spC1 " << spC1 << endl
      << " spC2 " << spC2 << endl
      << " iGauss " << iGauss << endl;

    if (method=="SN"){
        out_ << " bending_energy_ipart = Seung and Nelson" << endl;
        init_bendij(mesh);
        bending_energy_ipart = [this](Vec3d *pos, int *node_nbr, int num_nbr, int idx, int bdry_type, double lenth, int edge, bool pbc, double *lijsq) -> double {return this->SeungNelson(pos, node_nbr, num_nbr, idx, bdry_type,lenth, edge,pbc,lijsq);};

        out_ << "Bond based bending" << endl;
        exchange = [this](int idx1, int idx2, const MESH_p& mesh) -> void {
        return this->exchange_bond(idx1, idx2, mesh);};
    }else{
        out_ << " bending_energy_ipart = Itzykson" << endl;
        init_coefbend(mesh.compA, mesh.N);
        bending_energy_ipart = [this](Vec3d *pos, int *node_nbr, int num_nbr, int idx, int bdry_type, double lenth, int edge, bool pbc, double *lijsq) -> double {
        return this->Itzykson(pos, node_nbr, num_nbr, idx, bdry_type, lenth, edge,lijsq);};

        out_ << "Node based bending" << endl;
        exchange = [this](int idx1, int idx2, const MESH_p& mesh) -> void {
        return this->exchange_node(idx1, idx2);};
    }

    out_.close();
}
/*-------------------------------------------------*/
void BE::init_coefbend(int *lipA, int N){
    for (int i = 0; i < N; ++i) {
        if (lipA[i]) coef_bend.push_back(bend2);
        else coef_bend.push_back(bend1);
    }
}
/*--------------------------------------------------------------------------*/
void BE::init_bendij(MESH_p mesh) {
    int *lipA = mesh.compA;

    // Ensure `bendij` has the correct size based on `mesh.node_nbr_list`
    int totalNbrs = mesh.nghst * mesh.N;  // Or the maximum size of node_nbr_list
    bendij.resize(totalNbrs, 0.0);        // Resize instead of pushing

    for (int i = 0; i < mesh.N; ++i) {
        set_nbrbending(i, mesh);
        // double bendi = lipA[i] ? bend2 : bend1;  // Select bend2 if lipA[i] is true
        // int num_nbr = mesh.numnbr[i];
        // int cm_idx = mesh.nghst * i;

        // for (int k = cm_idx; k < cm_idx + num_nbr; ++k) {
        //     int j = mesh.node_nbr_list[k];

        //     // Validate `j` index to ensure it is within bounds of `lipA`
        //     if (j < 0 || j >= mesh.N) {
        //         std::cerr << "Error: Neighbor index " << j << " out of bounds." 
        //         << std::endl;
        //         continue;
        //     }

        //     double bendj = lipA[j] ? bend2 : bend1;  // Select bend2 for j if lipA[j] is true
        //     if (k < totalNbrs) {
        //         coef_bend[k] = (bendi + bendj) * 0.866;  // Store computed value
        //     } else {
        //         std::cerr << "Error: coef_bend index " << k << " out of bounds." << 
        //         std::endl;
        //     }
        // }
    }
}
/*------------------------------------------------------------------------------*/
void BE::set_nbrbending(int idx, const MESH_p& mesh){
    int *lipA = mesh.compA;
    int cm_idx = ghost*idx;
    // double bendi=coef_bend[idx];
    double bendi = lipA[idx] ? bend2 : bend1;
    int num_nbr = mesh.numnbr[idx];
    double bendj;
    for (int k = cm_idx; k < cm_idx+num_nbr; ++k) {
        int j = mesh.node_nbr_list[k];
        bendj = lipA[j] ? bend2 : bend1;
        bendij[k] = (bendi + bendj) * 0.866;  // Store computed value
    }
}
/*--------------------------------------------------------------------------*/
inline double acot(double x) {
    double result = atan(1.0 / x); // Calculate arccot(x)
    if (x < 0) {
        result += M_PI; // Adjust for the range [0, pi]
    }
    return result;
}
/*-------------------------------------------------*/
double BE::SeungNelson(Vec3d *pos, int *node_nbr, int num_nbr, int idx,
            int bdry_type, double lenth, int edge, bool pbc, double *lijsq){
    /// @brief Computes the bending energy contribution in Seung and Nelson way
    /// when the position of the ith particle changes.
    /// @param pos Array containing coordinates of all particles.
    /// @param idx Index of the ith particle.
    /// @param node_nbr Indices of nearest neighbors of the ith particle.
    /// @param num_nbr Number of nearest neighbors.
    /// @param bdry_type Type of boundary condition (unused).
    /// @param lijsq Preallocated array to store squared distances between neighbors.
    /// @return Bending energy contribution for the ith particle.
    /// @todo Merge calculations with SeungNelson_nbr for optimization using better
    /// data structures.
    /// @note Multiplies by 0.5 to avoid double-counting; 
    /// single-counted contributions cancel in energy differences.
    double bend_ener=0;
    Vec3d xij[num_nbr], ntri[num_nbr];
    int jdx, kdx;
    int nbrloopind;
    // this will decide whether we are circulating over all the neighbors or not. For boundary points (except pbc) we do not circulate over all the neighbors.
    if (idx>edge || pbc) {nbrloopind = num_nbr;}
    
    else {nbrloopind = num_nbr-1;}

    if (idx>edge || !pbc) for (int j = 0; j < nbrloopind; ++j){ xij[j]=pos[idx]-pos[node_nbr[j]]; }
    if (idx < edge && pbc) for (int j = 0; j < nbrloopind; ++j){ xij[j]=diff_pbc(pos[idx],pos[node_nbr[j]], lenth); }

    for (int j = 0; j < nbrloopind; ++j){lijsq[j]=normsq(xij[j]);}
    
    for (int j = 0; j < nbrloopind; ++j){ 
        ntri[j] = cross_product(xij[j],xij[(j+1)%num_nbr]);
        if (norm(ntri[j]) > 1e-10)  ntri[j] = ntri[j]/norm(ntri[j]);
    }
    for (int j = 0; j < nbrloopind; ++j){
        bend_ener+=bendij[idx*ghost+j]*(1-inner_product(ntri[j],ntri[(j+1)%num_nbr]));
    }
    return 0.5*bend_ener;
}
/*-------------------------------------------------*/
// There is a problem if we have the vertices on the polls -- zaxis.
// the code also needs to be fixed for pbc, boundary points etc.
// see SeungNelson function.
double BE::Itzykson(Vec3d *pos, int *node_nbr, int num_nbr, int idx,
            int bdry_type, double lenth, int edge, double *lijsq){
    /// @brief Estimate the Bending energy contribution when ith particle 
    /// position changes
    /// @param Pos array containing co-ordinates of all the particles
    /// @param idx index of ith particle;
    /// @param node_nbr nearest neigbours of idx; 
    /// @param num_nbr number of nearest neigbours of idx; 
    /// @param para  Membrane related parameters;
    /// @todo try openMP Pragmas;
    /// @return Bending Energy contribution when ith particle is displaced.
    double bend_ener, Gauss_ener=0;
    Vec3d cot_times_rij;
    Vec3d lap_bel,lap_bel_t0, nhat;
    Vec3d xjk;
    double cot_jdx_k,cot_kdx,cot_kmdx,area_ijkm;
    Vec3d xik,xikm,xjkm,nhat_local,xijp1;
    int jdx,kdx,kmdx,jdxp1;
    double liksq,likmsq,ljkmsq;
    double sigma_i = 0e0;
    double cot_aij[num_nbr],cot_bij[num_nbr],area_ijk[num_nbr];
    double ljksq[num_nbr];
    Vec3d xij[num_nbr];
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
        liksq=lijsq[(j+1)%num_nbr];
        likmsq=lijsq[(j-1+num_nbr)%num_nbr];
        area_ijkm=area_ijk[(j-1+num_nbr)%num_nbr];
        xijp1=xij[(j+1)%num_nbr];
        cot_times_rij = cot_times_rij + xij[j]*(cot_aij[j] + cot_bij[j]);
        sigma_i=sigma_i+voronoi_area(cot_jdx_k,cot_bij[j],liksq,lijsq[j],area_ijk[j]);
        nhat_local=cross_product(xijp1,xij[j]);
        nhat = nhat + nhat_local*(1e0/norm(nhat_local));
    }
    nhat = nhat/norm(nhat);
    lap_bel = cot_times_rij/sigma_i;
    lap_bel = lap_bel*0.5;
    lap_bel_t0 = nhat*spcurv;
    bend_ener = 0.5*sigma_i*normsq(lap_bel-lap_bel_t0);
    if (iGauss){
        Gauss_ener = M_PI*(2-num_nbr);
        for (int j = 0; j < num_nbr; ++j){
            Gauss_ener+=acot(cot_aij[j])+acot(cot_bij[j]);
        }
        bend_ener -= Gauss_ener;
    }
    // cout << bend_ener << endl;
    return coef_bend[idx]*bend_ener;
}
/*--------------------------------------------------------------------------*/
// Wrapper function
// double BE::bending_energy_ipart(Vec3d *pos, int *node_nbr, int num_nbr, int idx,
//             int bdry_type, double lenth, int edge){
//     double lijsq[num_nbr];
//     return bending_energy_ipart(pos, node_nbr, num_nbr, idx, bdry_type, lenth,
//         edge, lijsq);
// }
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
    double be=0e0;
    for(j = idx*mesh.nghst; j < idx*mesh.nghst + mesh.numnbr[idx]; j++){
        nbr = mesh.node_nbr_list[j];
        num_nbr_j = mesh.numnbr[nbr];
        cm_idx_nbr = nbr*mesh.nghst;
   	    double lijsq[num_nbr_j];
        be += bending_energy_ipart(pos,
            (int *) mesh.node_nbr_list + cm_idx_nbr,
            num_nbr_j, nbr, mesh.bdry_type, mesh.boxlen, mesh.lastbdry, mesh.pbc, lijsq);
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
    double be=0e0, lijsq[12];
    //
    st_idx = get_nstart(mesh.N, mesh.bdry_type);
    for(idx = st_idx; idx < mesh.N; idx++){
        /* idx = 2; */        
        cm_idx = idx*mesh.nghst;
        num_nbr = mesh.numnbr[idx];
        be+= bending_energy_ipart(pos, (int *) (mesh.node_nbr_list + cm_idx),
                num_nbr, idx, mesh.bdry_type, mesh.boxlen, mesh.lastbdry, mesh.pbc, lijsq);
    }
    return be;
}
/*------------------------------------------------------------------------------*/
void BE::exchange_node(int idx1, int idx2){
    swap(coef_bend[idx1], coef_bend[idx2]);
}
/*------------------------------------------------------------------------------*/
void BE::exchange_bond(int idx1, int idx2, const MESH_p& mesh){
    // It's important that you have already swiped the components.
    set_nbrbending(idx1, mesh);
    set_nbrbending(idx2, mesh);
}
/*----------------------------------------------------------------------------*/