#include <string>
#include <fstream>
#include "boundary.hpp"
#include "vector.hpp"
#include "mesh.hpp"

extern "C" void BdryRead(double*,bool*,bool*,bool*,bool*,bool*,char*);

int get_nstart(int,int);

int BDE::initBDE(int N, std::string fname){
  char tmp_fname[128];
  string parafile,outfile;

  parafile=fname+"/para_file.in";
  sprintf(tmp_fname,"%s",parafile.c_str());
  BdryRead(&coef_gc,&do_bdry,&fix_bottom,&fix_top,&fix_left,&fix_right,tmp_fname);

  ofstream out_;
  out_.open(fname+"/bdrypara.out");
  out_<< "#=========boundary parameters==========" << endl
      << " N " << N << endl
      << " coef_gc = " << coef_gc << endl
      << " do_bdry " << do_bdry << endl
      << " fix_bottom " << fix_bottom << endl
      << " fix_top " << fix_top << endl
      << " fix_left " << fix_left << endl
      << " fix_right " << fix_right << endl;
  out_.close();
  return 0;
}


//helper function definitions for geodesic curvature calculation
// NOTE: might be better to add in vector.cpp 

//finite difference derivative
//first derivative, three-point stencil
Vec3d dd1_3pt(Vec3d sm1, Vec3d s, Vec3d sp1){
  double h1=norm(s-sm1);
  double h2=norm(sp1-s);
  double denom=h1*h2*(h1+h2);
  Vec3d result;
  result.x=(-h2*sm1.x*h2 - (h1-h2)*s.x*(h1+h2) +h1*sp1.x*h1)/denom;
  result.y=(-h2*sm1.y*h2 - (h1-h2)*s.y*(h1+h2) +h1*sp1.y*h1)/denom;
  result.z=(-h2*sm1.z*h2 - (h1-h2)*s.z*(h1+h2) +h1*sp1.z*h1)/denom;  
  return result;
}

//second derivative, three-point stencil
Vec3d dd2_3pt(Vec3d sm1, Vec3d s, Vec3d sp1){
  double h1=norm(s-sm1);
  double h2=norm(sp1-s);
  double denom=h1*h2*(h1+h2);
  Vec3d result;
  result.x=2.0*(h2*sm1.x-(h1+h2)*s.x+h1*sp1.x)/denom;
  result.y=2.0*(h2*sm1.y-(h1+h2)*s.y+h1*sp1.y)/denom;
  result.z=2.0*(h2*sm1.z-(h1+h2)*s.z+h1*sp1.z)/denom;
  return result;
}

//get boundary vertices in positional order
void BDE::bdry_edges(Vec3d *pos, MESH_p mesh){
  if (!do_bdry) return;
  int nframe=get_nstart(mesh.N,mesh.bdry_type);
  //nframe is the number of boundary nodes
  num_bdry_nodes=nframe;
  
  //we had num_bdry_nodes=2*nghostx+2*nghosty (from init.cpp)
  //reconstruct nghostx and nghosty to isolate the edges

  n_bdry_x=0; //this is what was nghostx
  //we want to find the number of nodes in the bottom edge,
  //where the init.cpp placement started.
  //the (nghostx-1)th point, i.e. the (n_bdry_x-1)th point, is the first one in the nframe frame points having a neighbour list where i+1 should be missing.
  
  for (int i=0;i<num_bdry_nodes;i++){
    int cm=i*mesh.nghst;
    int nn=mesh.numnbr[i];
    int *nbrs=mesh.node_nbr_list+cm;
    bool next_found=false;
    for (int k=0; k<nn; k++){
      if (nbrs[k]==i+1){
        next_found=true;
        break;
      }
    }
    if (!next_found){
       n_bdry_x=i+1;
       break;
    }
  }
  n_bdry_y=(num_bdry_nodes-2*n_bdry_x)/2;

  //from init.cpp:
  //0...nghostx-1 was the bottom boundary
  //nghostx...2*nghostx-1 was the top boundary
  //2*nghostx...2*nghostx+nghosty-1 was the left boundary
  //2*nghostx+nghosty...2*nghostx+2*nghosty-1 was the right boundary
  
  is_clamped.assign(mesh.N,false);
  is_end.assign(mesh.N,false);
  which_edge.assign(mesh.N,-1);

  int edge_ends[]={n_bdry_x,2*n_bdry_x,2*n_bdry_x+n_bdry_y,num_bdry_nodes};
  int edge_starts[]={0,n_bdry_x,2*n_bdry_x,2*n_bdry_x+n_bdry_y};
  bool edge_fix_flag[]={fix_bottom,fix_top,fix_left,fix_right};

  for (int i=0;i<4;i++){ //go over each edge
    int start_idx=edge_starts[i];
    int end_idx=edge_ends[i]; 
 
    is_end[start_idx]=true;
    is_end[end_idx-1]=true; //mark the ends of edges
  
    //if we want to clamp one edge or something
    
    for (int j=start_idx;j<end_idx;j++){
        which_edge[j]=i;
        if (edge_fix_flag[i]){
           is_clamped[j]=true;
        }
    }
  }
}

//helper function, calculate normal for the edge
Vec3d BDE::midpt_normal(Vec3d *pos, MESH_p mesh,int v1, int v2){
  int cm1=v1*mesh.nghst;
  int nn1=mesh.numnbr[v1];
  int *nbrs1=mesh.node_nbr_list+cm1;

  int cm2=v2*mesh.nghst;
  int nn2=mesh.numnbr[v2];
  int *nbrs2=mesh.node_nbr_list+cm2;

  int v3=-1; //other vertex of the triangle, will be found
  for (int i=0; i<nn1; i++){
    for (int j=0; j<nn2; j++){
      if (nbrs1[i]==nbrs2[j]){
        v3=nbrs1[i]; //find the common neighbour (bdry edge has only one) 
        break;
      }
    }
    if (v3!=-1) break;
  }
  
  if (v3==-1){ //if error - third pt not found
    Vec3d defaultnormal={0.0,0.0,1.0}; //direct outward z
    return defaultnormal;
  }  

  Vec3d edgeA=pos[v2]-pos[v1];
  Vec3d edgeB=pos[v3]-pos[v1];
  Vec3d N=cross_product(edgeA,edgeB);
   
  //check for normal orientation
  //for flat ribbon normal direction should be outwards (up?)
  //this needs to be changed for twisting etc.
  if (N.z < 0) {
    N.x = -N.x;
    N.y = -N.y;
    N.z = -N.z;
  }
 	 
  double Nnorm=norm(N);
  Vec3d N_unit;
  N_unit.x=N.x/Nnorm;
  N_unit.y=N.y/Nnorm;
  N_unit.z=N.z/Nnorm;
  return N_unit;
}

 
//skip endpoints for now  

double BDE::midpt_kg(Vec3d *pos, MESH_p mesh, int v1, int v2){
  if (is_end[v1]||is_end[v2]) return 0.0;
  
  Vec3d pvert1=pos[v1-1];
  Vec3d vert1=pos[v1];
  Vec3d nvert1=pos[v1+1];

  Vec3d pvert2=pos[v2-1];
  Vec3d vert2=pos[v2];
  Vec3d nvert2=pos[v2+1];

  Vec3d T1=dd1_3pt(pvert1,vert1,nvert1);
  Vec3d dTds1=dd2_3pt(pvert1,vert1,nvert1);

  Vec3d T2=dd1_3pt(pvert2,vert2,nvert2);
  Vec3d dTds2=dd2_3pt(pvert2,vert2,nvert2);

  //we are assuming v1,v2 are not endpoints, i.e. the 3pt stencils exist.
  //we are getting the tangent and dtds at v1 and v2
  //now we interpolate at the midpoint

  Vec3d T=Vec3d_add(T1,T2,1.0);
  Vec3d dTds=Vec3d_add(dTds1,dTds2,1.0);

  dTds.x=0.5*dTds.x;
  dTds.y=0.5*dTds.y;
  dTds.z=0.5*dTds.z;

  //T will be normalised anyways so we can skip the 0.5* multiplication

  double Tnorm=norm(T);
  Vec3d T_unit;
  T_unit.x=T.x/Tnorm;
  T_unit.y=T.y/Tnorm;
  T_unit.z=T.z/Tnorm;

  //we have unit tangent and dTds 
  Vec3d N_unit=midpt_normal(pos,mesh,v1,v2);
  
  Vec3d NcrossT=cross_product(N_unit,T_unit);
  double kg=inner_product(dTds,NcrossT);
  //this is the geodesic curvature at the midpoint

  return kg;
}

double BDE::bde_ipart(Vec3d *pos, MESH_p mesh, int idx){
  if (!do_bdry||idx>=num_bdry_nodes) return 0.0;
  double E=0.0; //we will find the energy contribution

  //affected edges when idx is moved 
  int  edge_list[4][2]={{idx-2,idx-1},{idx-1,idx},{idx,idx+1},{idx+1,idx+2}};

  for (int i=0; i<4;i++){
    int v1=edge_list[i][0];
    int v2=edge_list[i][1];

    if (v1>=0 && v2<num_bdry_nodes && !is_end[v1] && !is_end[v2]){
      double kg=midpt_kg(pos,mesh,v1,v2);
      double ds=norm(pos[v2]-pos[v1]);
      E-=coef_gc*kg*ds;
      //coef_gc is the saddle-splay modulus
      //the energy contribution of the Gaussian curvature is the NE>
      //we subtract the contribution due to kg
      //the Euler characteristic is treated as an irrelevant topological constant 

    }
  }
  return E;
}

double BDE::bde_total(Vec3d *pos, MESH_p mesh){
  if (!do_bdry) return 0.0;
  double E=0.0;
  for (int idx=0; idx<num_bdry_nodes; idx++){
    if (is_end[idx]) continue;
    if (idx-1 >=0 && !is_end[idx-1]){
      double kg=midpt_kg(pos,mesh,idx-1,idx);
      double ds=norm(pos[idx]-pos[idx-1]);
      E-=coef_gc*kg*ds;
    }
  }
  return E;
}


bool BDE::is_clamped_vertex(int idx){ //to be passed into metropolis
  if (idx>=num_bdry_nodes) return false;
  return is_clamped[idx];
}

void BDE::initBdryCache(Vec3d *pos, MESH_p mesh){
  if (!do_bdry) return;
  bdry_cache.assign(mesh.N,0.0);
  for (int idx=0;idx<num_bdry_nodes;idx++){
    bdry_cache[idx]=bde_ipart(pos,mesh,idx);
  }
}
