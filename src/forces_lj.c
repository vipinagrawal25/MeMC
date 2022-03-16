#include "../include/global.h"
#include "../include/subroutine.h"
 
double cal_length(double x1 , double x2, 
        double y1, double y2, double len, char
        *metric){
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

    /* if (metric == "cart" ){ */
    /* printf("cartesian %s %d %d \n", metric, is_sph, is_cart); */
    /* } */
    if (is_cart == 1 ){
        dx = x2 - x1;
        dy = y2 - y1;
        if(dx >  len*0.5) dx = dx - len; 
        if(dy >  len*0.5) dy = dy - len;
        if(dx < -len*0.5) dx = dx + len;
        if(dy < -len*0.5) dy = dy + len;
        ds = dx*dx + dy*dy;
    }
    if (is_sph == 1 ){
        // x = theta
        // y is the phi
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


void make_nlist_pf(POSITION *Pos, Neighbours *neib,
        LJpara para, char *metric){

    int cellx,celly;
    double r2_cut, new_rc;
    double displ; 
    int cnt;


    r2_cut = para.r_cut*para.r_cut;

    for(int i=0;i < para.N; i++)
        neib[i].cnt_ss = 0;


    displ = 0.0;
    cnt = 0;

    double cellcut = para.r_cut + R_del; 
    new_rc = cellcut*cellcut;

    for(int i = 0; i < para.N; i++){
        for(int j = i+1; j < para.N; j++){
            int m = i;
            int n = j;
            if(len_check_pf(Pos[m], Pos[n],  para.len, new_rc, metric)){
                neib[m].list_ss[neib[m].cnt_ss] = n;
                neib[n].list_ss[neib[n].cnt_ss] = m;
                neib[m].cnt_ss += 1;
                neib[n].cnt_ss += 1;
            }
        }

    }
}


bool len_check_pf(POSITION s1, POSITION s2,
        double len, double new_rc, char *metric){

    double ds;
    bool var = false;

    ds = cal_length(s1.x , s2.x, 
            s1.y, s2.y, len, metric);
    var = (ds < new_rc);
    return var;
}

bool len_check_ss(POSITION s1, POSITION s2, 
        double len, double new_rc){

    double drx = s1.x - s2.x;
    double dry = s1.y - s2.y;
    bool var = false;

    if(drx > len*0.5) drx = drx - len;
    if(dry > len*0.5) dry = dry - len;
    if(drx < - len*0.5) drx = drx + len;
    if(dry < - len*0.5) dry = dry + len;
    double Sq_dr = drx*drx + dry*dry;
    var = (Sq_dr < new_rc);
    return var;
}

double pairlj_ipart_energy(POSITION *Pos, int *n_list,
        int ni, int i_p, LJpara para, char *metric){
    double Sq_dr2;
    int j, k;
    double dx,dy,Sq_dr1,r2,r6;
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
                loc_ener += 4.0*eps*r6*(r6);
            }
        }
    }
    return loc_ener;
}


double pairlj_total_energy(POSITION *Pos, Neighbours *neib,
        LJpara para, char *metric){
    int i,j,k,ni;
    double j_pe, pe;

    pe = 0e0;
    for(i=0;i<para.N;i++){
        ni = neib[i].cnt_ss;
        j_pe = pairlj_ipart_energy(Pos, neib[i].list_ss,
                ni, i, para, metric);
        pe += j_pe;
    }
    return pe;
}
