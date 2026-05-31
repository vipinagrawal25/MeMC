#include "subroutine.h"
#include "global.h"
using namespace std;
void init_eval_lij_t0(Vec3d *Pos, MESH mesh, double *lij_t0, int N);
void init_mbrane(MBRANE_para *mbrane);
//
int main(int argc, char const *argv[]){
	string folder=argv[1];
	//
	string prefix="part_";
	int fnum=1;
	string extension;
	Vec3d *Pos;
	MESH mesh;
    ofstream bend_fid(folder+"/bendE.txt", ios::out);
  	ofstream stretch_fid(folder+"/stretchE.txt", ios::out);
	double *lij_t0;
	MBRANE_para mbrane;
	int num_nbr, cm_idx, idx;
	double be,se;
	init_mbrane(&mbrane);
	// allocate arrays
    Pos = (Vec3d *)calloc(mbrane.N, sizeof(Vec3d));
	mesh.cmlist = (int *)calloc(mbrane.N+1, sizeof(int));
    mesh.node_nbr_list = (int *)calloc(mbrane.num_nbr, sizeof(int));
    lij_t0 = (double *)calloc(mbrane.num_nbr, sizeof(double));
    //
	hdf5_io_read_pos( (double *)Pos, folder+"/input.h5");
	hdf5_io_read_mesh((int *) mesh.cmlist,
	        (int *) mesh.node_nbr_list, folder+"/input.h5");
	init_eval_lij_t0(Pos, mesh, lij_t0, mbrane.N);
	//
	if(FileExists(folder+"/"+prefix+ZeroPadNumber(fnum)+".vtk"))extension=".vtk";
	else if(FileExists(folder+"/"+prefix+ZeroPadNumber(fnum)+".h5"))extension=".h5";
	string fname=folder+"/"+prefix+ZeroPadNumber(fnum)+extension;
	//
	while(FileExists(fname)){
		if (extension==".vtk"){
			visit_vtk_read((double*) Pos, mbrane.N, fname);
            // print((double*)Pos, 3*mbrane.N);
		}else if(extension==".h5"){
			hdf5_io_read_pos((double*) Pos, fname);
		}
		for(idx = 0; idx < mbrane.N; idx++){
	        num_nbr = mesh.cmlist[idx + 1] - mesh.cmlist[idx];
	        cm_idx = mesh.cmlist[idx];
	        be  =bending_energy_ipart(Pos, 
	                (int *) (mesh.node_nbr_list + cm_idx),
	                 num_nbr, idx, mbrane);
	        se  =stretch_energy_ipart(Pos, 
                (int *) (mesh.node_nbr_list + cm_idx),
                (double *) (lij_t0 + cm_idx), num_nbr, 
                idx, mbrane );
	        bend_fid << be << "\n";
	        stretch_fid << se << "\n";
    	}
    	bend_fid << "\n";
    	stretch_fid << "\n";
        fnum=fnum+1;
        fname=folder+"/"+prefix+ZeroPadNumber(fnum)+extension;
        cout << fname << endl;
    }
    bend_fid.close();
    stretch_fid.close();
	return 0;
}
void init_eval_lij_t0(Vec3d *Pos, MESH mesh, double *lij_t0, int N){
    /// @brief evaluates distance between neighbouring points and stores in lij_t0
    ///  @param Pos array containing co-ordinates of all the particles
   /// @param lij_t0 initial distance between points of membrane
    ///  @param mesh mesh related parameters -- connections and neighbours
    /// information; 
    ///  @param para membrane related parameters 
    Vec3d dr;
    int i,j,k;
    int num_nbr, cm_idx;
    double sum_lij=0;
    double r0;
    for(i = 0; i < N; i++){
        num_nbr = mesh.cmlist[i + 1] - mesh.cmlist[i];
        cm_idx = mesh.cmlist[i];
        for(k = cm_idx; k < cm_idx + num_nbr; k++) {
            j = mesh.node_nbr_list[k];
            dr.x = Pos[i].x - Pos[j].x;
            dr.y = Pos[i].y - Pos[j].y;
            dr.z = Pos[i].z - Pos[j].z;
            lij_t0[k] = sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
            sum_lij+=lij_t0[k];
        }
    }
}
void init_mbrane(MBRANE_para *mbrane){
	mbrane->N=5120;
	mbrane->coef_bend=4096;
	mbrane->YY=4616*mbrane->coef_bend;       // Young's modulus
	//
    mbrane->coef_vol_expansion=0;
    mbrane->sp_curv = 2;
    mbrane->pressure = 0;
    mbrane->radius = 1;
    mbrane->pos_bot_wall = 0;
    mbrane->sigma = 0;
    mbrane->epsilon = 0;
    mbrane->theta = 0;
    mbrane->num_triangles = 2*mbrane->N - 4;
    mbrane->num_nbr = 3*mbrane->num_triangles;
}