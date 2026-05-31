#include "global.h"
#include "subroutine.h"

extern "C" void Membrane_listread(int *, double *, double *, 
        double *, int *, char *);

extern "C" void Spcurv_listread(char*, double *, double *, double *, char *);

extern "C" void Stick_listread(bool *, double *, double *, double *, 
        double *,  char *);

extern "C" void  MC_listread(char *, double *, double *, bool *,
           int *, int *, char *);

extern "C"  void  Activity_listread(char *, double *, double *, char *);

extern "C"  void  Afm_listread(bool *, double *, double *, double *, double *,
             char *);

extern "C"  void  Spring_listread(bool *, int *, double *, double *, char *);

extern "C"  void  Fluid_listread(bool *, int * , int *, double *, char *);
extern "C" void   Volume_listread(bool *, bool *, double *, double*, char *); 
extern "C" void   Area_listread(bool *, double *, char *); 

/*--------------------------------------------------------------------------------*/
void init_eval_lij_t0(Vec3d *Pos, MESH_p mesh, double *lij_t0,
         MBRANE_p *para, SPRING_p *spring, bool is_fluid){
    /// @brief evaluates distance between neighbouring points and stores in lij_t0
    ///  @param Pos array containing co-ordinates of all the particles
   /// @param lij_t0 initial distance between points of membrane
    ///  @param mesh mesh related parameters -- connections and neighbours
    /// information; 
    ///  @param para membrane related parameters 
    Vec3d dr;
    int i,j,k;
    int num_nbr, cm_idx, npairs;
    double sum_lij=0;
    double r0;
    npairs = 0;
    for(i = 0; i < para->N; i++){
        num_nbr = mesh.numnbr[i];
        cm_idx = mesh.nghst * i;
        for(k = cm_idx; k < cm_idx + num_nbr; k++) {
            j = mesh.node_nbr_list[k];
           /* dr = diff_pbc(Pos[i] , Pos[j], para->len); */
            dr = Pos[j] - Pos[i];
            lij_t0[k] = sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
            sum_lij += sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
            npairs++;
            /* printf("%g %g %g %g %g \n", Pos[i].x, Pos[j].x, Pos[i].y, Pos[j].y, lij_t0[k]); */
        }
    }
    para->av_bond_len = sum_lij/npairs;
    r0=para->av_bond_len;
    spring->constant=para->coef_bend/(r0*r0);
    if(is_fluid){
        for(i = 0; i < mesh.nghst*para->N; i++){
            lij_t0[i] = para->av_bond_len;
        }
    }
}
/*--------------------------------------------------------------------------*/
void init_area_t0(Vec3d *pos, MESH_p mesh, MBRANE_p mbrane_para, AREA_p area_para){
    /// @brief Compute area for each triangle
    // int Nt = 2*mbrane_para.N-4;
    int num_nbr,cm_idx, jdx, jdxp1;
    Vec3d xij, xijp1;
    double area_tot=0;
    int st_idx = get_nstart(mbrane_para.N, mbrane_para.bdry_type);
    for(int idx = st_idx; idx < mbrane_para.N; idx++){
        /* idx = 2; */
        num_nbr = mesh.numnbr[idx];
        cm_idx = mesh.nghst*idx;
        area_ipart(pos, (double *) (area_para.area_t0 + cm_idx),
                  (int *) (mesh.node_nbr_list + cm_idx),
                  num_nbr, idx);
        for (int k = cm_idx+num_nbr; k < cm_idx+mesh.nghst; ++k){
            area_para.area_t0[k] = -1;
        }
    }
    // print(area_para.area_t0,mesh.nghst*mbrane_para.N);
}
/*--------------------------------------------------------------------------*/
bool check_param(MBRANE_p *mbrane_para, SPCURV_p *spcurv_para, MC_p *mc_para, 
        FLUID_p *fld_para, VOL_p *vol_para, AREA_p *area_para, STICK_p *stick_para, 
        AFM_p *afm_para,  ACTIVE_p *act_para, SPRING_p *spring_para){
    ///@brief This function checks all the parameters and warns the user 
    /// in case of unphysical parameters.
    bool status=true;
    if (fabs(vol_para->coef_vol_expansion)<1e-8){vol_para->do_volume=false;}
    if (fabs(vol_para->pressure)<1e-8){vol_para->is_pressurized=false;}
    if (fabs(area_para->coef_area_expansion)<1e-8){area_para->do_area=false;}
    // if (!(mc_para->algo=="mpolis")||!(mc_para->algo=="glauber")){
    //     cout<< "I do not understand your choice of mc algo." << endl
    //     << "setting it to mpolis" << endl;
    //     mc_para->algo="mpolis";
    // }
    if (vol_para->do_volume==true && vol_para->is_pressurized==true){
        cout<< "The shell can not be pressured while conserving the volume." << endl
        << "Set one of do_volume or is_pressurized equal to zero" << endl;
        status=false;
    }
    return status;
}
/*--------------------------------------------------------------------------*/
bool init_read_parameters(MBRANE_p *mbrane_para, SPCURV_p *spcurv_para, MC_p *mc_para, 
        FLUID_p *fld_para, VOL_p *vol_para, AREA_p *area_para,
        STICK_p *stick_para, AFM_p *afm_para,  ACTIVE_p *act_para, 
        SPRING_p *spring_para, string para_file){
   /// @brief read parameters from para_file 
    ///  @param mesh mesh related parameters -- connections and neighbours
    /// information; 
    ///  @param mbrane membrane related parameters
    ///  @param afm AFM related parameters
    ///  @param mcpara monte-carlo related parameters
    //
    //
    char temp_algo[char_len];
    char which_act[char_len], tmp_fname[char_len];
    char spcurv_which[char_len];
    //
    sprintf(tmp_fname, "%s", para_file.c_str() );
    Membrane_listread(&mbrane_para->N, &mbrane_para->coef_bend,
            &mbrane_para->YY, &mbrane_para->radius,
            &mbrane_para->bdry_type, tmp_fname);

    // sprintf(tmp_fname, "%s", para_file.c_str() );
    sprintf(tmp_fname, "%s", para_file.c_str() );
    Spcurv_listread(spcurv_which, &spcurv_para->minC, &spcurv_para->maxC,
            &spcurv_para->theta, tmp_fname);
    spcurv_para->which_spcurv=spcurv_which;

    sprintf(tmp_fname, "%s", para_file.c_str() );
    Stick_listread(&stick_para->do_stick, &stick_para->pos_bot_wall, 
            &stick_para->sigma, &stick_para->epsilon, &stick_para->theta,
            tmp_fname);

    sprintf(tmp_fname, "%s", para_file.c_str());
    MC_listread(temp_algo, &mc_para->dfac, &mc_para->kBT, &mc_para->is_restart,
            &mc_para->tot_mc_iter, &mc_para->dump_skip, tmp_fname);
    mc_para->algo=temp_algo;

    sprintf(tmp_fname, "%s", para_file.c_str() );
    Spring_listread(&spring_para->do_spring, &spring_para->icompute, &spring_para->nPole_eq_z,
            &spring_para->sPole_eq_z, tmp_fname);

    sprintf(tmp_fname, "%s", para_file.c_str() );
    Afm_listread(&afm_para->do_afm, &afm_para->tip_rad, 
            &afm_para->tip_pos_z, &afm_para->sigma, &afm_para->epsilon,
             tmp_fname);

    sprintf(tmp_fname, "%s", para_file.c_str() );
    Activity_listread(which_act, &act_para->minA, &act_para->maxA, tmp_fname);
    act_para->act = which_act;

    // sprintf(tmp_fname, "%s", para_file.c_str() );
    sprintf(tmp_fname, "%s", para_file.c_str() );
    Fluid_listread(&fld_para->is_fluid, &fld_para->min_allowed_nbr,
            &fld_para->fluidize_every, &fld_para->fac_len_vertices, tmp_fname);

    sprintf(tmp_fname, "%s", para_file.c_str() );
    Volume_listread(&vol_para->do_volume, &vol_para->is_pressurized,
            &vol_para->coef_vol_expansion, &vol_para->pressure, tmp_fname);

    sprintf(tmp_fname, "%s", para_file.c_str());
    Area_listread(&area_para->do_area, &area_para->coef_area_expansion, tmp_fname);
    // mbrane->av_bond_len = sqrt(8*pi/(2*mbrane->N-4));
    // define the monte carlo parameters
    mc_para->one_mc_iter = 2*mbrane_para->N;
    mc_para->delta = sqrt(8*pi/(2*mbrane_para->N-4));
    // string filename = "00000/para.out";
    return true;
}
//
void write_parameters(MBRANE_p mbrane, SPCURV_p spcurv_para, MC_p mc_para,
        FLUID_p fld_para, VOL_p vol_p, AREA_p area_p, STICK_p stick_para, AFM_p afm_para,  
        ACTIVE_p act_para, SPRING_p spring_para, string out_file){
 
    double FvK = mbrane.YY*mbrane.radius*mbrane.radius/mbrane.coef_bend;
    ofstream out_;
    out_.open( out_file );
    out_<< "# =========== Model Parameters ==========" << endl
            << " Foppl von Karman number: FvK = " << FvK << endl
            << " Elasto-thermal number: ET = " << mc_para.kBT/mbrane.coef_bend*sqrt(FvK) << endl
            << " average bond length: r0 = " << mbrane.av_bond_len << endl;

    out_<< "# =========== Membrane Parameters ==========" << endl
            << " coef bend = " << mbrane.coef_bend << endl
            << " YY " << mbrane.YY << endl
            // << " sp_curv " << mbrane.sp_curv << endl
            << " N " << mbrane.N << endl
            << " av_bond_len " << mbrane.av_bond_len << endl
            << " bdry_type " << mbrane.bdry_type << endl
            << " radius " << mbrane.radius << endl;

    out_<< "# =========== Spontaneous Curvature Parameters ==========" << endl
            << " which Curvature = " << spcurv_para.which_spcurv << endl
            << " minC = " << spcurv_para.minC << endl
            << " maxC = " << spcurv_para.maxC << endl
            << " theta =" << spcurv_para.theta << endl;

    out_<< "# =========== Monte Carlo Parameters ==========" << endl
            << " algo = " << mc_para.algo << endl
            << " dfac " << mc_para.dfac << endl
            << " dump_skip " << mc_para.dump_skip << endl
            << " kbt " << mc_para.kBT << endl
            << " delta " << mc_para.delta << endl
            << " is_restart " << mc_para.is_restart << endl
            << " tot_mc_iter " << mc_para.tot_mc_iter << endl
            << " one mc iter " << mc_para.one_mc_iter << endl;


    out_<< "# =========== Activity Parameters ==========" << endl
            << " which activity = " << act_para.act << endl
            << " minA = " << act_para.minA << endl
            << " maxA = " << act_para.maxA << endl;


    out_<< "# =========== Fluid Parameters ==========" << endl
            << " is fluid= " << fld_para.is_fluid << endl
            << " min_allowed_nbr = " << fld_para.min_allowed_nbr << endl
            << " fluid iter every " << fld_para.fluidize_every << endl
            << " factor_len_vertices = " << fld_para.fac_len_vertices << endl;

    out_<< "# =========== Volume Parameters ==========" << endl
            << " do volume= " << vol_p.do_volume << endl
            << " is pressurized = " << vol_p.is_pressurized << endl
            << " coef_vol_expansion " << vol_p.coef_vol_expansion << endl
            << " pressure  " << vol_p.pressure << endl;

    out_<< "# =========== Area Parameters ==========" << endl
            << " do area = " << area_p.do_area << endl
            << " coef_area_expansion = " << area_p.coef_area_expansion << endl;

    out_<< "# =========== Sticking Parameters ==========" << endl
            << " do stick " << stick_para.do_stick << endl
            << " pos_bot_wall  " << stick_para.pos_bot_wall << endl
            << " sigma " << stick_para.sigma << endl
            << " epsilon " << stick_para.epsilon << endl
            << " theta " << stick_para.theta << endl;

    out_<< "# =========== Spring Parameters ==========" << endl
            << " do spring " << spring_para.do_spring << endl
            << " icompute  " << spring_para.icompute << endl
            << " constant " << spring_para.constant << endl
            << " nPole_eq_z " << spring_para.nPole_eq_z << endl
            << " sPole_eq_z " << spring_para.sPole_eq_z << endl;

    out_.close();
}
/*--------------------------------------------------------------------------*/
void init_activity(ACTIVE_p activity, int N){
    int i;

    if(activity.act == "random"){
        for(i=0;i<N;i++) activity.activity[i] = RandomGenerator::generateUniform(activity.minA, activity.maxA);
    }
    if(activity.act == "constant"){
        for(i=0;i<N;i++) activity.activity[i] = activity.maxA;
    }
}
/*--------------------------------------------------------------------------*/
void init_spcurv(SPCURV_p spcurv, Vec3d *pos, int N){
    int i;
    double theta;
    if(spcurv.which_spcurv=="delta"){
        spcurv.spcurv[0]=spcurv.maxC;
        for (int i = 1; i < N; ++i){spcurv.spcurv[i]=spcurv.minC;}
    }
    int temp=0;
    if(spcurv.which_spcurv=="constant"){
        for(i= 0; i<N; i++){
            theta = pi - acos(pos[i].z);
            temp+= (theta<spcurv.theta);
            if (theta<spcurv.theta){spcurv.spcurv[i]=spcurv.maxC;}
            else{spcurv.spcurv[i]=spcurv.minC;}
        }
    }
}
/*****************************************************************************/
