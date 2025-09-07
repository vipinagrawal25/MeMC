/**
 * @file flat.hpp
 * @brief Preprocessing the cells and vertices before it passes to the code.
 *        So that the code runs efficiently.
 *        Step: 1) Make sure that the boundary points are at the boundary
 *              2) From pos and cells -> pos, cumu_list, node_nbr. Make sure to handle both formats.
 *              3) The neighbours and boundary points should be stored in clockwise manner.
 *              4) The lastbdry returns the last points of the boundary even for periodic mesh.
 *                  a) Basically aperiodic and periodic meshes would return the same number in that case.
 */
#ifndef FLAT_HPP
#define FLAT_HPP

#include <vector>
#include <algorithm>
#include <set>
#include <utility>
#include "vector.hpp"
#include <iostream>
#include <cmath>
#include <cassert>
#include "misc.hpp"

set<pair<int, int>> make_bond_list(int *node_nbr, int N, int nghst, int *numnbr);
using namespace std;

void get_neighbours(int *node_nbr, int *numnbr, const vector<vector<int>> &simplices, int N, int nghst){
    // Initialize all entries in node_nbr to -1
    fill(node_nbr, node_nbr + N * nghst, -1);

    // Initialize numnbr array to 0
    fill(numnbr, numnbr + N, 0);

    // For each simplex (triangle), add neighbors for each vertex
    for (const auto &simplex : simplices){
        // For each vertex in the triangle, add the other two vertices as neighbors
        for (int i = 0; i < 3; ++i){
            int v = simplex[i];            // Current vertex
            int v1 = simplex[(i + 1) % 3]; // Next vertex
            int v2 = simplex[(i + 2) % 3]; // Previous vertex

            // Check if v1 is already in the neighbor list of v
            bool v1_found = false;
            bool v2_found = false;
            for (int j = 0; j < numnbr[v]; ++j){
                if (node_nbr[v * nghst + j] == v1)
                    v1_found = true;
                if (node_nbr[v * nghst + j] == v2)
                    v2_found = true;
            }

            // Add v1 if not already present and within limit
            if (!v1_found && numnbr[v] < nghst){
                node_nbr[v * nghst + numnbr[v]] = v1;
                numnbr[v]++;
            }

            // Add v2 if not already present and within limit
            if (!v2_found && numnbr[v] < nghst){
                node_nbr[v * nghst + numnbr[v]] = v2;
                numnbr[v]++;
            }
        }
    }
}

vector<vector<int>> sort_simplices(int *cells, int num_cells){
    // Equivalent to Python's sort_simplices function
    vector<vector<int>> sorted_cells;
    sorted_cells.reserve(num_cells * 6); // Pre-allocate space for efficiency
    // For each cell, create all 6 possible orientations
    for (int i = 0; i < num_cells; i++)
    {
        vector<int> scles = {cells[i * 3], cells[i * 3 + 1], cells[i * 3 + 2]};
        sort(scles.begin(), scles.end()); // Sort the original triangle vertices
        // Add all 6 possible orientations
        sorted_cells.push_back({scles[0], scles[1], scles[2]});
        sorted_cells.push_back({scles[1], scles[2], scles[0]});
        sorted_cells.push_back({scles[2], scles[0], scles[1]});
        sorted_cells.push_back({scles[0], scles[2], scles[1]});
        sorted_cells.push_back({scles[1], scles[0], scles[2]});
        sorted_cells.push_back({scles[2], scles[1], scles[0]});
    }

    sort(sorted_cells.begin(), sorted_cells.end(),
         [](const vector<int> &a, const vector<int> &b)
         {
             return (a[0] < b[0]) ||
                    (a[0] == b[0] && a[1] < b[1]) ||
                    (a[0] == b[0] && a[1] == b[1] && a[2] < b[2]);
         });
    return sorted_cells;
}

vector<vector<int>> remove_duplicates(vector<vector<int>> &sorted_cells){
    auto last = unique(sorted_cells.begin(), sorted_cells.end());
    sorted_cells.erase(last, sorted_cells.end());
    return sorted_cells;
}

set<pair<int, int>> make_bond_list(int *node_nbr, int N, int nghst, int *numnbr){
    set<pair<int, int>> edge_set;
    // Determine the number of vertices from the length of neighbor_indices and nghst.
    for (int i = 0; i < N; ++i){
        int current_index = i * nghst;
        // Loop over the valid neighbors for vertex i
        for (int j = 0; j < numnbr[i]; ++j){
            int neighbor = node_nbr[current_index + j];
            // Create the edge (i, neighbor) in sorted order
            pair<int, int> edge = make_pair(min(i, neighbor), max(i, neighbor));
            edge_set.insert(edge);
        }
    }
    return edge_set;
}
/*-------------------*/
bool is_wrapped_edge(int n1, int n2, const Vec3d* pos, double boxL) {
    double dx = pos[n1].x - pos[n2].x;
    double dy = pos[n1].y - pos[n2].y;
    double dz = pos[n1].z - pos[n2].z;

    // Apply minimum image convention
    double dx_mi = dx - boxL * round(dx / boxL);
    double dy_mi = dy - boxL * round(dy / boxL);
    double dz_mi = dz - boxL * round(dz / boxL);

    // If original != min-image, it’s a wrapped edge
    double d2_orig = dx*dx + dy*dy + dz*dz;
    double d2_mi   = dx_mi*dx_mi + dy_mi*dy_mi + dz_mi*dz_mi;

    return fabs(d2_orig - d2_mi) > 1e-8;
}

bool is_boundary_vertex(Vec3d *pos, const vector<int> &nbrs, const set<pair<int, int>> &edge_set, double length){
    /// @brief The function checks if a vertex is a boundary vertex by examining its neighbors and the edge set. We do not care whether it's PBC or not. Both returns the same result.
    /// That is, if a vertex is a boundary vertex in a PBC mesh, it will also be a boundary vertex in an aperiodic mesh and vice versa.
    /// @param nbrs List of neighbors for the vertex
    /// @param edge_set Set of edges in the mesh
    /// @return true if the vertex is a boundary vertex, false otherwise

    size_t n = nbrs.size();
    if (n < 3){
        return true;
    }
    for (size_t i = 0; i < n; ++i){
        int n1 = nbrs[i];
        int n2 = nbrs[(i + 1) % n];
        // Create an edge with sorted order
        pair<int, int> edge = make_pair(min(n1, n2), max(n1, n2));
        // If the edge is not found in the edge set, it's a boundary edge
        if (edge_set.find(edge) == edge_set.end()){
            return true;
        }

        // Edge exists but is periodic-wrapped ⇒ ignore, treat as missing
        if (is_wrapped_edge(n1, n2, pos, length)){
            return true;
        }

    }
    return false;
}

// Sort an index list of points in CCW order around their centroid (uses x,y of Vec3d)
inline void sort_indices_ccw_by_centroid(vector<int> &indices, const Vec3d *pos){
    if (indices.empty()) return;
    double cx = 0.0, cy = 0.0;
    for (int idx : indices){
        cx += pos[idx].x;
        cy += pos[idx].y;
    }
    cx /= static_cast<double>(indices.size());
    cy /= static_cast<double>(indices.size());

    vector<pair<double,int>> ang_idx;
    ang_idx.reserve(indices.size());
    for (int idx : indices){
        double dx = pos[idx].x - cx;
        double dy = pos[idx].y - cy;
        double ang = atan2(dy, dx);
        ang_idx.emplace_back(ang, idx);
    }
    sort(ang_idx.begin(), ang_idx.end(), [](const pair<double,int> &a, const pair<double,int> &b){
        return a.first < b.first;
    });
    // write back
    for (size_t i = 0; i < ang_idx.size(); ++i) indices[i] = ang_idx[i].second;

    // rotate so that the ordering starts at a deterministic boundary point
    // choose the point with minimum x (tie-breaker min y) as the anchor
    int anchor_pos = 0;
    for (size_t i = 1; i < indices.size(); ++i){
        int cur = indices[i];
        int best = indices[anchor_pos];
        if (pos[cur].x < pos[best].x || (pos[cur].x == pos[best].x && pos[cur].y < pos[best].y)){
            anchor_pos = static_cast<int>(i);
        }
    }
    if (anchor_pos != 0){
        rotate(indices.begin(), indices.begin() + anchor_pos, indices.end());
    }
}

bool is_all_points_on_one_side(Vec3d v1, Vec3d v2, const vector<Vec3d> &points, 
    bool verbose = false){
    // Check if all neighbor points are on one side of the line formed by the first two
    Vec3d line_vec = v2 - v1;
    line_vec = line_vec/norm(line_vec);
    vector<double> cross_zs;
    for(auto point : points){
        auto vec = point - v1;
        auto cross_z = line_vec.x * vec.y - line_vec.y * vec.x;
        cross_zs.push_back(cross_z);
    }
    if(verbose){
        cout << "Cross products: ";
        for (auto z : cross_zs) cout << z << " ";
        cout << endl;
    }
    bool all_positive = all_of(cross_zs.begin(), cross_zs.end(), [](double z){ return z >= 0; });
    bool all_negative = all_of(cross_zs.begin(), cross_zs.end(), [](double z){ return z <= 0; });
    return all_positive || all_negative;
}

inline vector<Vec3d> vec_of_nbr_points(Vec3d *pos, const vector<int> &nbrs, size_t start_idx, size_t end_idx){
    vector<Vec3d> nbr_points;
    for (size_t k = start_idx; k < end_idx; ++k){
        nbr_points.push_back(pos[nbrs[k]]);
    }
    return nbr_points;
}

void order_boundary_neighbors(int *node_nbr, int *numnbr, Vec3d *pos, int Nb, int nghst = 12){
    // Rearrange the boundary neighbors in a consistent order
    // such that it starts at the boundary and ends at the boundary.
    // This is important for correctly
    // Nb: number of boundary vertices
    // pos: array of vertex positions (Vec3d), indexed by global vertex id
    // numnbr: array containing number of neighbors for each vertex
    // node_nbr: flat neighbor list (length >= Nb*nghst), neighbors for vertex i start at node_nbr[nghst*i]
    bool logic;
    for (int ib = 0; ib < Nb; ++ib){
        int num_nbr = numnbr[ib];
        if (num_nbr <= 0)
            continue;
        int start = ib * nghst;
        // vector so you can use it in the rotate function.
        vector<int> nbrs;
        nbrs.reserve(num_nbr);
        for (int k = 0; k < num_nbr; ++k)
            nbrs.push_back(node_nbr[start + k]);
        //
        logic = is_all_points_on_one_side(pos[ib], pos[nbrs[0]], vec_of_nbr_points(pos, nbrs, 1, num_nbr)) && is_all_points_on_one_side(pos[ib], pos[nbrs[num_nbr-1]], vec_of_nbr_points(pos, nbrs, 0, num_nbr-1));
        for (int attempt = 0; !logic && attempt < num_nbr; ++attempt){
            // rotate right by one element
            rotate(nbrs.begin(), nbrs.end() - 1, nbrs.end());
            logic = is_all_points_on_one_side(pos[ib], pos[nbrs[0]], vec_of_nbr_points(pos, nbrs, 1, num_nbr)) && is_all_points_on_one_side(pos[ib], pos[nbrs[num_nbr - 1]], vec_of_nbr_points(pos, nbrs, 0, num_nbr - 1));
        }
        if (!logic){
            cerr << "Error: Could not order neighbors of boundary vertex " << ib << " at (" << pos[ib].x << ", " << pos[ib].y << ", " << pos[ib].z << ") with " << num_nbr << " neighbors." << endl;
            exit(EXIT_FAILURE);
        }
    }
}

int put_boundary_first(Vec3d *pos, int *node_nbr_list, int *numnbr, int N, int nghst, double length){
    /// @brief This function makes sure that the boundary points are at the beginning of the list. Steps: 1) make bond list 2) identify boundary points 3) rearrange the points such that boundary points are at the beginning. 4) rearrange the neighbours accordingly.
    /// @param pos Array of vertex positions
    /// @param numnbr Array containing the number of neighbours for each vertex
    /// @param node_nbr_list Array containing the list of neighbours for each vertex
    /// @param N Total number of vertices
    /// @param nghst Maximum number of neighbours per vertex
    /// @return The last index of boundary vertices

    auto allbonds = make_bond_list(node_nbr_list, N, nghst, numnbr);

    vector<int> boundary_points;
    vector<int> interior_points;
    
    for (int i = 0; i < N; ++i){
        int *nbrs = node_nbr_list + i * nghst;
        vector<int> nbr_list(nbrs, nbrs + numnbr[i]);
        if (is_boundary_vertex(pos, nbr_list, allbonds, length)){
            boundary_points.push_back(i);
            // printf("Boundary vertex %d at (%.3f, %.3f, %.3f) with %d neighbors\n", i, pos[i].x, pos[i].y, pos[i].z, numnbr[i]);
        }
        else{
            interior_points.push_back(i);
        }
    }
    // print(boundary_points);
    // Now rearrange pos, numnbr, node_nbr_list such that boundary points are at the beginning.
    // Sort boundary points in CCW order around their centroid to provide a consistent ordering.
    if (!boundary_points.empty()){
        sort_indices_ccw_by_centroid(boundary_points, pos);
    }

    vector<int> new_order;
    new_order.insert(new_order.end(), boundary_points.begin(), boundary_points.end());
    new_order.insert(new_order.end(), interior_points.begin(), interior_points.end());
    
    // Create mapping from old indices to new indices
    vector<int> old_to_new(N);
    for (int i = 0; i < N; ++i){
        old_to_new[new_order[i]] = i;
    }
    //
    Vec3d *new_pos = new Vec3d[N];
    int *new_numnbr = new int[N];
    int *new_node_nbr_list = new int[N * nghst];
    //
    for (int i = 0; i < N; ++i){
        int old_index = new_order[i];
        new_pos[i] = pos[old_index];
        new_numnbr[i] = numnbr[old_index];
        // Rearrange the neighbour list and update indices
        for (int j = 0; j < nghst; ++j){
            int old_nbr = node_nbr_list[old_index * nghst + j];
            if (old_nbr >= 0 && old_nbr < N) {
                new_node_nbr_list[i * nghst + j] = old_to_new[old_nbr];
            } else {
                new_node_nbr_list[i * nghst + j] = old_nbr; // Keep invalid indices as is
            }
        }
    }
    // // Copy back to original arrays
    copy(new_pos, new_pos + N, pos);
    copy(new_numnbr, new_numnbr + N, numnbr);
    copy(new_node_nbr_list, new_node_nbr_list + N * nghst, node_nbr_list);
    //
    delete[] new_pos;
    delete[] new_numnbr;
    delete[] new_node_nbr_list;
    // 
    return boundary_points.empty() ? -1 : static_cast<int>(boundary_points.size()) - 1;
}

/**
 * @brief Sort neighbors in flat (2D) configuration with periodic boundary conditions
 * @param pts Array of 2D points (assumed to have x, y coordinates)
 * @param Np Number of points
 * @param numnbr Array containing number of neighbors for each point
 * @param node_nbr Input/output neighbor array 
 * @param nghst Maximum neighbors per point
 * @param length Periodic box length for boundary condition handling
 */
void sort_nbrs(int* node_nbr, Vec3d* pts, int Np, int* numnbr, int nghst, double length){
    for (int ip = 0; ip < Np; ++ip){
        int num_nbr = numnbr[ip];
        int start_idx = ip * nghst;
        
        if (num_nbr <= 1) continue; // No need to sort if 0 or 1 neighbors
        
        // Store angles and corresponding neighbor indices
        vector<pair<double, int>> angle_nbr_pairs;
        
        for (int j = 0; j < num_nbr; ++j) {
            int nbr_idx = node_nbr[start_idx + j];
            
            // Calculate dx and dy with periodic boundary conditions
            double dx = pts[nbr_idx].x - pts[ip].x;
            double dy = pts[nbr_idx].y - pts[ip].y;
            
            // Handle periodic boundary conditions
            if (dx >= 0.5 * length) dx -= length;
            if (dx < -0.5 * length) dx += length;
            if (dy >= 0.5 * length) dy -= length;
            if (dy < -0.5 * length) dy += length;
            
            // Calculate angle
            double angle = atan2(dy, dx);
            
            // Convert negative angles to positive (0 to 2π range)
            if (angle < 0) {
                angle += 2.0 * M_PI;
            }    
            angle_nbr_pairs.push_back({angle, nbr_idx});
        }
        // Sort by angle
        sort(angle_nbr_pairs.begin(), angle_nbr_pairs.end());
        // Update node_nbr with sorted neighbors
        for (int j = 0; j < num_nbr; ++j) {
            node_nbr[start_idx + j] = angle_nbr_pairs[j].second;
        }
    }
}

// // Helper: check if all points (given by indices) lie on one side of the line from v1 to v2
// inline bool is_all_points_on_one_side(const Vec3d &v1, const Vec3d &v2, const int *point_indices, int num_points, const Vec3d *pts){
//     const double eps = 1e-12; // tolerance for numerical precision
//     // line vector in 2D
//     double lx = v2.x - v1.x;
//     double ly = v2.y - v1.y;
//     double norm = sqrt(lx*lx + ly*ly);
//     if (norm <= eps) return true; // degenerate line — treat as "one side"
//     lx /= norm; ly /= norm;

//     // compute signed cross products for each point
//     bool all_nonneg = true;
//     bool all_nonpos = true;
//     for (int i = 0; i < num_points; ++i){
//         const Vec3d &p = pts[point_indices[i]];
//         double vx = p.x - v1.x;
//         double vy = p.y - v1.y;
//         // 2D cross product (line_vec x vec)
//         double cross = lx * vy - ly * vx;
//         if (cross < -eps) all_nonneg = false;
//         if (cross > eps) all_nonpos = false;
//         if (!all_nonneg && !all_nonpos) return false;
//     }
//     return all_nonneg || all_nonpos;
// }

// // Reorders boundary neighbors so that they are oriented/rolled into a consistent ordering.
// // Renamed from `sort_bdry` to `order_boundary_neighbors` for clarity.
// // Nb    : number of boundary vertices
// // pts   : array of vertex positions (Vec3d), indexed by global vertex id
// // cmlist: cumulative list (int array) describing neighbor ranges per boundary vertex
// // node_nbr: flat neighbor list (length >= Nb*nghst), neighbors for vertex i start at node_nbr[nghst*i]
// // nghst : stride (max neighbors per vertex)
// inline void order_boundary_neighbors(int Nb, const Vec3d *pts, const int *numnbr, int *node_nbr, int nghst = 12){
//     for (int ib = 0; ib < Nb; ++ib){
//         // number of neighbors for this boundary vertex
//         int num_nbr = numnbr[ib];
//         if (num_nbr <= 0) continue;

//         int start = nghst * ib;
//         // gather neighbor indices
//         vector<int> nbrs;
//         nbrs.reserve(num_nbr);
//         for (int k = 0; k < num_nbr; ++k) nbrs.push_back(node_nbr[start + k]);
        
//         auto check_logic = [&](const vector<int> &arr)->bool{
//             if (arr.size() < 2) return true;
//             // first check: pts[ib] -> pts[arr[0]] against arr[1:]
//             if (!is_all_points_on_one_side(pts[ib], pts[arr[0]], arr.data()+1, static_cast<int>(arr.size()-1), pts)) return false;
//             // second check: pts[ib] -> pts[arr.back()] against arr[0:-1]
//             if (!is_all_points_on_one_side(pts[ib], pts[arr.back()], arr.data(), static_cast<int>(arr.size()-1), pts)) return false;
//             return true;
//         };
//         bool logic = check_logic(nbrs);
//         // try rolling the neighbour list until condition satisfied or until we've tried all rotations
//         for (int attempt = 0; !logic && attempt < num_nbr; ++attempt){
//             // rotate right by one element
//             rotate(nbrs.begin(), nbrs.end()-1, nbrs.end());
//             logic = check_logic(nbrs);
//         }
//         if (!logic) {
//             cerr << "order_boundary_neighbors: failed to order neighbors for boundary vertex " << ib
//                       << " — keeping original order.\n";
//             exit(EXIT_FAILURE);
//         }
//         // write back the (possibly rotated) neighbours
//         for (int k = 0; k < num_nbr; ++k) node_nbr[start + k] = nbrs[k];
//     }
// }

/**
 * Python function for sorting neighbors in anticlockwise order
 * 
 * def sort_nbrs(R, Np, cmlst, node_nbr):
 *     zhat = np.array([0.,0.,1.])
 *     for i in range(Np):
 *         nbrs=node_nbr[cmlst[i]:cmlst[i+1]]  # neighbours of ith node
 *         vector=R[i]
 *         # I will rotate the coordinate system about this vector
 *         vhat = np.cross(vector,zhat)       
 *         vnorm = LA.norm(vhat)
 *         # If the vector is already lying at z-axis then there is no need to rotate.
 *         if vnorm>1e-16:
 *             vhat = vhat/vnorm
 *             theta = polar(vector)
 *             # Rotate all the neighbours of a point.
 *             rotated=rotate(R[nbrs],vhat,theta)
 *             # Since all the voronoi cells are rotated, sort them in anticlockwise direction
 *             sorted_indices = sort_2Dpoints_theta(rotated[:,0],rotated[:,1])[0]
 *             node_nbr[cmlst[i]:cmlst[i+1]]=nbrs[sorted_indices]
 *     return node_nbr
 */

#endif