#include "stretching.hpp"
#include <fstream>
#include "misc.hpp"

extern "C" void StretchRead(double *, bool *, bool *, 
                            double *, double *, double *, bool *, char *);

int get_nstart(int, int);

STE::STE(const MESH_p& mesh, std::string fname){
    char tmp_fname[128];
    double radius=mesh.radius;
    string parafile, outfile;

    parafile = fname+"/para_file.in";
    sprintf(tmp_fname, "%s", parafile.c_str() );
    StretchRead(&YY, &do_volume, &is_pressurized, &Kappa,
              &pressure, &coef_area_expansion, &do_area, tmp_fname);
    ini_vol = mesh.ini_vol;
    cout << "ini_vol=" << ini_vol << endl;
    cout << do_volume << endl;
    ofstream out_;
    out_.open( fname+"/stretchpara.out");
    out_<< "# =========== stretching parameters ==========" << endl
      << " N " << mesh.N << endl
      << " coef_stretch = " << YY*sqrt(3)/2 << endl
      << " do_volume " << do_volume << endl
      << " is_pressurized " << is_pressurized << endl
      << " pressure " << pressure << endl
      << " coef_vol_expansion " << Kappa << endl
      << " do_area " << do_area << endl
      << " coef_area_expansion " << coef_area_expansion << endl;
    out_.close();
}

inline Vec3d diff(Vec3d a, Vec3d b, double lenth, int bdry_type, int idx, int edge){
    if (bdry_type==1||idx<=edge) return a - b;
    else    return diff_pbc(a, b, lenth);
}

double STE::stretch_energy_ipart(Vec3d *pos,int *node_nbr, int num_nbr, int idx, 
            int ghost, int bdry_type, double lenth, int edge){
   double HH;
   double idx_ener;
   Vec3d rij;
   double mod_rij, avg_lij;
   int i,j;
   //
   idx_ener = 0e0;
   // HH = para.coef_str/(para.av_bond_len*para.av_bond_len);
   HH = YY*sqrt(3)/2;
   if (bdry_type == 1 || idx>edge){
      for (i =0; i < num_nbr; i++){
         j = node_nbr[i];
         rij = pos[idx] - pos[j];
         // rij = diff(pos[idx], pos[j], lenth, bdry_type, idx, edge);
         mod_rij = sqrt(inner_product(rij, rij));
         idx_ener = idx_ener + (mod_rij - lij_t0[idx*ghost+i])*(mod_rij - lij_t0[idx*ghost+i]);
      }
   }else{
      for (i =0; i < num_nbr; i++){
         j = node_nbr[i];
         rij = diff_pbc(pos[idx], pos[j], lenth);
         mod_rij = sqrt(inner_product(rij, rij));
         idx_ener = idx_ener + (mod_rij - lij_t0[idx*ghost+i])*(mod_rij - lij_t0[idx*ghost+i]);
      }
   }
   return 0.5*idx_ener*HH;
}

double STE::stretch_energy_total(Vec3d *pos, MESH_p mesh){
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
    st_idx = get_nstart(mesh.N, mesh.bdry_type);
    se = 0e0;
    for(idx = st_idx; idx < mesh.N; idx++){
        /* idx = 2; */
        num_nbr = mesh.numnbr[idx];
        cm_idx = idx*mesh.nghst;
        se += stretch_energy_ipart(pos, (int *)(mesh.node_nbr_list + cm_idx),
                num_nbr, idx, mesh.nghst, mesh.bdry_type, mesh.boxlen, mesh.edge);
    }
    return se*0.5e0;
}
/*----------------------------------------------------------------------*/
double STE::area_ipart(Vec3d *pos, int *node_nbr, int num_nbr, int idx, 
            int bdry_type, double lenth, int edge){
    ///@brief This function computes the
    ///area of all the triangles near the vertices.
    int jdx, jdxp1;
    Vec3d xij, xijp1;
    double area_energy_idx=0;
    double area_idx=0;
    if (bdry_type==1||idx>edge){
      for (int k = 0; k < num_nbr; ++k){
        jdx = node_nbr[k];
        jdxp1 = node_nbr[(k+1)%num_nbr];
        xij = pos[idx]- pos[jdx];
        xijp1 = pos[idx] - pos[jdxp1];
        area_idx = area_idx + 0.5*norm(cross_product(xij,xijp1));
      }
    }else{
      for (int k = 0; k < num_nbr; ++k){
        jdx = node_nbr[k];
        jdxp1 = node_nbr[(k+1)%num_nbr];
        xij = change_vec_pbc(pos[idx]-pos[jdx], lenth);
        xijp1 = change_vec_pbc(pos[idx]-pos[jdxp1],lenth);
        area_idx = area_idx + 0.5*norm(cross_product(xij,xijp1));
      }
    }
    return area_idx;
}
/*----------------------------------------------------------------------*/
double STE::area_total(Vec3d *pos, MESH_p mesh){
     int idx, st_idx;
     int num_nbr, cm_idx;
     double area;

     st_idx = get_nstart(mesh.N, mesh.bdry_type);
     area = 0e0;
     for(idx = st_idx; idx < mesh.N; idx++){
         /* idx = 2; */
         cm_idx = idx*mesh.nghst;
         num_nbr = mesh.numnbr[idx];
         area += area_ipart(pos,  
                 (int *) (mesh.node_nbr_list + cm_idx),
                 num_nbr, idx, mesh.bdry_type, mesh.boxlen, mesh.edge);
     }
     return area/3;
}
/*----------------------------------------------------------------------*/
double STE::init_eval_lij_t0(MESH_p &mesh, bool is_fluid){
    /// @brief evaluates distance between neighbouring points and stores in lij_t0
    ///  @param Pos array containing co-ordinates of all the particles
   ///  @param lij_t0 initial distance between points of membrane
    ///  @param mesh mesh related parameters -- connections and neighbours
    /// information; 
    ///  @param para membrane related parameters 
    Vec3d dr;
    int i,j,k;
    int num_nbr, cm_idx, npairs;
    double sum_lij=0;
    double tlij;
    double r0, av_bond_len;
    npairs = 0;
    Vec3d *Pos = mesh.pos;
    double lenth = mesh.boxlen;

    for(i = 0; i < mesh.nghst*mesh.N; i++){
        lij_t0.push_back(0.0);
    }

    for(i = 0; i < mesh.N; i++){
        num_nbr = mesh.numnbr[i];
        cm_idx = mesh.nghst * i;
        for(k = cm_idx; k < cm_idx + num_nbr; k++) {
            j = mesh.node_nbr_list[k];
            dr = diff_pbc(Pos[j], Pos[i], lenth);
            lij_t0[k] = sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
            sum_lij += sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
            npairs++;
            /* printf("%g %g %g %g %g \n", Pos[i].x, Pos[j].x, Pos[i].y, Pos[j].y, lij_t0[k]); */
        }
    }
    // print(mesh.node_nbr_list,mesh.N*mesh.nghst);
    av_bond_len = sum_lij/npairs;
    // spring->constant=para->coef_bend/(r0*r0);
    if(is_fluid){
        for(i = 0; i < mesh.nghst*mesh.N; i++){
            lij_t0[i] = av_bond_len;
        }
    }
    // print(mesh.node_nbr_list,mesh.N*mesh.nghst);
    return av_bond_len;
}

double STE::volume_ipart(Vec3d *pos, int *node_nbr,
        int num_nbr, int idx, int bdry_type, double lenth, int edge){
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
    //
    volume1 = 0e0;
    ri = pos[idx];
    if (bdry_type==1||idx>edge){
      for (i =0; i < num_nbr; i++){
        j = node_nbr[i];
        k=node_nbr[(i+1)%num_nbr];
        rj = pos[j]; rk = pos[k];
        rijk = (ri + rj + rk)*1/3e0;
        rij  = ri -  rj;
        rjk  = ri - rk;
        area1 = cross_product(rjk, rij);
        double ip1 = 0.5*inner_product(area1,rijk);
        volume1 = volume1 + abs(ip1);
      }
    }else{
      for (i =0; i < num_nbr; i++){
        j = node_nbr[i];
        k=node_nbr[(i+1)%num_nbr];
        rj = pos[j]; rk = pos[k];
        rijk=change_vec_pbc(ri+rj+rk, lenth);
        rij=change_vec_pbc(ri -  rj, lenth);
        rjk  = change_vec_pbc(ri -  rk, lenth);
        area1 = cross_product(rjk, rij);
        double ip1 = 0.5*inner_product(area1,rijk);
        volume1 = volume1 + abs(ip1);
      }
    }
    volume1 = volume1/3e0;
    return volume1;
}


double STE::volume_total(Vec3d *pos, MESH_p mesh){
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
     for(idx = 0; idx < mesh.N; idx++){
         /* idx = 2; */
        cm_idx = idx*mesh.nghst;
        num_nbr = mesh.numnbr[idx];
        vol += volume_ipart(pos,
                 (int *) (mesh.node_nbr_list + cm_idx),
                  num_nbr, idx, mesh.bdry_type, mesh.boxlen, mesh.edge);
     }
     return vol/3e0;
}

double STE::vol_energy_change(double vi, double dvol){
    double vf = vi+dvol;
    double de_vol = (vf/ini_vol - 1)*(vf/ini_vol - 1) - (vi/ini_vol - 1)*(vi/ini_vol - 1);
    // double de_vol = (2*dvol/(ini_vol*ini_vol))*(volume  - ini_vol) 
    //                 + (dvol/ini_vol)*(dvol/ini_vol);
    return Kappa*de_vol;
}
// /*-------------------------------------------------------------------------------------*/
// double area_energy_ipart(Vec3d *pos, int *node_nbr, double *area_t0, int num_nbr,
//                         int idx, double coef_area_expansion){
//     int jdx, jdxp1;
//     Vec3d xij, xijp1;
//     double area_energy_idx=0;
//     double area[num_nbr];
//     area_ipart(pos, area, node_nbr, num_nbr, idx);
//     for (int k = 0; k < num_nbr; ++k){
//         area_energy_idx += (1 - area[k]/area_t0[k])*(1 - area[k]/area_t0[k]);
//     }
//     return coef_area_expansion/2*(area_energy_idx);
// }
/*-------------------------------------------------------------------------------------*/

// double area_energy_total(Vec3d *pos, MESH_p mesh, MBRANE_p para, AREA_p area_para){
//     int idx, st_idx;
//     int num_nbr, cm_idx;
//     double ae;
//     st_idx = get_nstart(para.N, para.bdry_type);
//     ae = 0e0;
//     for(idx = st_idx; idx < para.N; idx++){
//         num_nbr = mesh.numnbr[idx];
//         cm_idx = idx*mesh.nghst;
//         ae += area_energy_ipart(pos,
//                 (int *)(mesh.node_nbr_list + cm_idx),
//                 (double *) (area_para.area_t0 + cm_idx),
//                 num_nbr, idx, area_para.coef_area_expansion);
//     }
//     return ae/3e0;
// }
/*-------------------------------------------------------------------------------------*/
