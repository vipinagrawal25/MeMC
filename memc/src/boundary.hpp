#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include <string>
#include <vector>
#include "mesh.hpp"
#include "vector.hpp"

class BDE{
  public:
    int initBDE(int N, std::string fname);
    void bdry_edges(Vec3d *pos, MESH_p mesh);
    void initBdryCache(Vec3d *pos, MESH_p mesh);

    Vec3d midpt_normal(Vec3d *pos, MESH_p mesh, int v1, int v2);
    double midpt_kg(Vec3d *pos, MESH_p mesh, int v1, int v2);
    double bde_ipart(Vec3d *pos, MESH_p mesh, int idx);
    double bde_total(Vec3d *pos, MESH_p mesh);

    bool is_clamped_vertex(int idx);

    bool do_bdry;
    double coef_gc;
    int num_bdry_nodes;
    int n_bdry_x;
    int n_bdry_y;

    std::vector<bool> is_clamped;
    std::vector<bool> is_end;
    std::vector<int> which_edge;
    std::vector<int> bdry_cache;

  private:
    bool fix_bottom;
    bool fix_top;
    bool fix_left;
    bool fix_right;
};
Vec3d dd1_3pt(Vec3d sm1, Vec3d s, Vec3d sp1);
Vec3d dd2_3pt(Vec3d sm1, Vec3d s, Vec3d sp1);
#endif
