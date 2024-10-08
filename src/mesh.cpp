#include "mesh.hpp"
using namespace std;
void MESH::initMESH(){
    MeshRead(&bdry_type, &nghst, &radius, tmp_fname);
    N = (int)hdf5_io_get_Np(outfolder+"/input.h5", "pos")/3;
    Pos = (Vec3d *)calloc(N, sizeof(Vec3d));
    numnbr = (int *)calloc(N, sizeof(int));
    node_nbr_list = (int *)calloc(nghst*N, sizeof(int));
}