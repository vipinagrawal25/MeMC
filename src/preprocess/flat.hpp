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
#include <map>
#include <utility>
#include "vector.hpp"
#include "mesh.hpp"
#include <iostream>
#include <cmath>
#include <cassert>
#include "misc.hpp"

set<pair<int, int>> make_bond_list(int *node_nbr, int N, int nghst, int *numnbr);
using namespace std;
vector<vector<int>> sort_simplices(int *cells, int num_cells){
    // Equivalent to Python's sort_simplices function
    vector<vector<int>> sorted_cells;
    sorted_cells.reserve(num_cells * 6); // Pre-allocate space for efficiency
    // For each cell, create all 6 possible orientations
    for (int i = 0; i < num_cells; i++){
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
    
    // Apply repeat_unique_pairs functionality
    // Create pairs from first two columns and count occurrences
    map<pair<int, int>, int> pair_counts;
    
    for (const auto &row : sorted_cells) {
        if (row.size() >= 2) {
            pair<int, int> p = make_pair(row[0], row[1]);
            pair_counts[p]++;
        }
    }
    
    // Find unique pairs (count == 1)
    set<pair<int, int>> unique_pairs_set;
    for (const auto &entry : pair_counts) {
        if (entry.second == 1) {
            unique_pairs_set.insert(entry.first);
        }
    }
    
    // Create new result vector with both original and repeated pairs in sorted order
    vector<vector<int>> result;
    result.reserve(sorted_cells.size() + unique_pairs_set.size());
    
    for (const auto &row : sorted_cells) {
        // Add the original row
        result.push_back(row);
        
        // Check if this row's pair (first two elements) is a unique pair
        if (row.size() >= 2) {
            pair<int, int> current_pair = make_pair(row[0], row[1]);
            if (unique_pairs_set.count(current_pair) > 0) {
                // Insert the repeated pair with -1 right after this row
                vector<int> repeated_row = {row[0], row[1], -1};
                result.push_back(repeated_row);
            }
        }
    }
    
    return result;
}

vector<vector<int>> remove_duplicates(vector<vector<int>> &sorted_cells){
    auto last = unique(sorted_cells.begin(), sorted_cells.end());
    sorted_cells.erase(last, sorted_cells.end());
    return sorted_cells;
}
/*-------------------*/
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

// bool is_boundary_vertex(Vec3d *pos, const vector<int> &nbrs, const set<pair<int, int>> &edge_set, double length){
//     /// @brief The function checks if a vertex is a boundary vertex by examining its neighbors and the edge set. We do not care whether it's PBC or not. Both returns the same result.
//     /// That is, if a vertex is a boundary vertex in a PBC mesh, it will also be a boundary vertex in an aperiodic mesh and vice versa.
//     /// @param nbrs List of neighbors for the vertex
//     /// @param edge_set Set of edges in the mesh
//     /// @return true if the vertex is a boundary vertex, false otherwise

//     size_t n = nbrs.size();
//     if (n < 3){
//         return true;
//     }
//     for (size_t i = 0; i < n; ++i){
//         int n1 = nbrs[i];
//         int n2 = nbrs[(i + 1) % n];
//         // Create an edge with sorted order
//         pair<int, int> edge = make_pair(min(n1, n2), max(n1, n2));
//         // If the edge is not found in the edge set, it's a boundary edge
//         if (edge_set.find(edge) == edge_set.end()){
//             return true;
//         }

//         // Edge exists but is periodic-wrapped ⇒ ignore, treat as missing
//         if (is_wrapped_edge(n1, n2, pos, length)){
//             return true;
//         }

//     }
//     return false;
// }

// Sort an index list of points in CCW order around their centroid (uses x,y of Vec3d)
// inline void sort_indices_ccw_by_centroid(vector<int> &indices, const Vec3d *pos){
//     if (indices.empty()) return;
//     double cx = 0.0, cy = 0.0;
//     for (int idx : indices){
//         cx += pos[idx].x;
//         cy += pos[idx].y;
//     }
//     cx /= static_cast<double>(indices.size());
//     cy /= static_cast<double>(indices.size());

//     vector<pair<double,int>> ang_idx;
//     ang_idx.reserve(indices.size());
//     for (int idx : indices){
//         double dx = pos[idx].x - cx;
//         double dy = pos[idx].y - cy;
//         double ang = atan2(dy, dx);
//         ang_idx.emplace_back(ang, idx);
//     }
//     sort(ang_idx.begin(), ang_idx.end(), [](const pair<double,int> &a, const pair<double,int> &b){
//         return a.first < b.first;
//     });
//     // write back
//     for (size_t i = 0; i < ang_idx.size(); ++i) indices[i] = ang_idx[i].second;

//     // rotate so that the ordering starts at a deterministic boundary point
//     // choose the point with minimum x (tie-breaker min y) as the anchor
//     int anchor_pos = 0;
//     for (size_t i = 1; i < indices.size(); ++i){
//         int cur = indices[i];
//         int best = indices[anchor_pos];
//         if (pos[cur].x < pos[best].x || (pos[cur].x == pos[best].x && pos[cur].y < pos[best].y)){
//             anchor_pos = static_cast<int>(i);
//         }
//     }
//     if (anchor_pos != 0){
//         rotate(indices.begin(), indices.begin() + anchor_pos, indices.end());
//     }
// }

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
//
bool is_boundary_vertex(const vector<int> &nbrs, const set<pair<int, int>> &edge_set){
    size_t n = nbrs.size();
    if (n < 3){
        return true;
    }
    for (size_t i = 0; i < n; ++i){
        int n1 = nbrs[i];
        int n2 = nbrs[(i + 1) % n];
        pair<int, int> edge = make_pair(min(n1, n2), max(n1, n2));
        if (edge_set.find(edge) == edge_set.end()){
            return true;
        }
    }
    return false;
}

// Overload that takes positions and box length so periodic (wrapped) edges
// are also treated as boundary edges when appropriate. If pos==nullptr or
// boxL<=0 the behavior is identical to the two-argument overload.
bool is_boundary_vertex(const vector<int> &nbrs, const set<pair<int, int>> &edge_set,
                        const Vec3d *pos, double boxL){
    size_t n = nbrs.size();
    if (n < 3){
        return true;
    }
    for (size_t i = 0; i < n; ++i){
        int n1 = nbrs[i];
        int n2 = nbrs[(i + 1) % n];
        pair<int, int> edge = make_pair(min(n1, n2), max(n1, n2));
        // missing edge -> boundary
        if (edge_set.find(edge) == edge_set.end()){
            return true;
        }
        // if geometry provided, also check whether the edge is periodic-wrapped
        // (minimum-image differs from original). If so, treat as boundary here.
        if (pos != nullptr && boxL > 0.0){
            if (is_wrapped_edge(n1, n2, pos, boxL)){
                // cout << "Wrapped edge detected between " << n1 << " and " << n2 << endl;
                return true;
            }
        }
    }
    return false;
}
//
void order_boundary_neighbors(int *node_nbr, int *numnbr, Vec3d *pos, 
    const set<pair<int, int>> &edge_set, int N, int nghst = 12){
    // Rearrange the boundary neighbors in a consistent order
    // such that it starts at the boundary and ends at the boundary.
    // This is important for correctly
    // Nb: number of boundary vertices
    // pos: array of vertex positions (Vec3d), indexed by global vertex id
    // numnbr: array containing number of neighbors for each vertex
    // node_nbr: flat neighbor list (length >= Nb*nghst), neighbors for vertex i start at node_nbr[nghst*i]
    bool logic;
    for (int i = 0; i < N; ++i){
        int num_nbr = numnbr[i];
        if (num_nbr <= 0) continue;
        int start = i * nghst;
        int attempt;
        int *nbradd = node_nbr + start;
        vector<int> nbrs(nbradd, nbradd + num_nbr);
        if (is_boundary_vertex(nbrs, edge_set)){
            logic = is_all_points_on_one_side(pos[i], pos[nbrs[0]], vec_of_nbr_points(pos, nbrs, 1, num_nbr)) && is_all_points_on_one_side(pos[i], pos[nbrs[num_nbr-1]], vec_of_nbr_points(pos, nbrs, 0, num_nbr-1));
            for (attempt = 0; !logic && attempt < num_nbr; ++attempt) {
                rotate(nbrs.begin(), nbrs.end() - 1, nbrs.end());
                logic = is_all_points_on_one_side(pos[i], pos[nbrs[0]], vec_of_nbr_points(pos, nbrs, 1, num_nbr)) && is_all_points_on_one_side(pos[i], pos[nbrs[num_nbr - 1]], vec_of_nbr_points(pos, nbrs, 0, num_nbr - 1));
            }
            for (int k = 0; k < num_nbr; ++k)
                node_nbr[start + k] = nbrs[k];
            if (!logic){
                cerr << "Error: Could not order neighbors of boundary vertex " << i << " at (" << pos[i].x << ", " << pos[i].y << ", " << pos[i].z << ") with " << num_nbr << " neighbors." << endl;
                exit(EXIT_FAILURE);
            }
        }
    }
}
//
set<pair<int, int>> make_bond_list(int *node_nbr, int *numnbr, int N, int nghst) {
    set<pair<int, int>> edge_set;
    // Determine the number of vertices from the length of neighbor_indices and nghst.
    for (int i = 0; i < N; ++i) {
        int current_index = i * nghst;
        // Loop over the valid neighbors for vertex i
        for (int j = 0; j < numnbr[i]; ++j){
            int neighbor = node_nbr[current_index + j];
            pair<int, int> edge = make_pair(min(i, neighbor), max(i, neighbor));
            edge_set.insert(edge);
        }
    }
    return edge_set;
}
// /**
//  * @brief Sort neighbors in flat (2D) configuration with periodic boundary conditions
//  * @param pts Array of 2D points (assumed to have x, y coordinates)
//  * @param Np Number of points
//  * @param numnbr Array containing number of neighbors for each point
//  * @param node_nbr Input/output neighbor array 
//  * @param nghst Maximum neighbors per point
//  * @param length Periodic box length for boundary condition handling
//  */
// void sort_nbrs(int* node_nbr, Vec3d* pts, int Np, int* numnbr, int nghst, double length) {
//     for (int ip = 0; ip < Np; ++ip){
//         int num_nbr = numnbr[ip];
//         int start_idx = ip * nghst;
        
//         if (num_nbr <= 1) continue; // No need to sort if 0 or 1 neighbors
        
//         // Store angles and corresponding neighbor indices
//         vector<pair<double, int>> angle_nbr_pairs;
        
//         for (int j = 0; j < num_nbr; ++j) {
//             int nbr_idx = node_nbr[start_idx + j];
            
//             // Calculate dx and dy with periodic boundary conditions
//             double dx = pts[nbr_idx].x - pts[ip].x;
//             double dy = pts[nbr_idx].y - pts[ip].y;
            
//             // Handle periodic boundary conditions
//             if (dx >= 0.5 * length) dx -= length;
//             if (dx < -0.5 * length) dx += length;
//             if (dy >= 0.5 * length) dy -= length;
//             if (dy < -0.5 * length) dy += length;
            
//             // Calculate angle
//             double angle = atan2(dy, dx);
//             // Convert negative angles to positive (0 to 2π range)
//             if (angle < 0) {
//                 angle += 2.0 * M_PI;
//             }    
//             angle_nbr_pairs.push_back({angle, nbr_idx});
//         }
//         // Sort by angle
//         sort(angle_nbr_pairs.begin(), angle_nbr_pairs.end());
//         // Update node_nbr with sorted neighbors
//         for (int j = 0; j < num_nbr; ++j) {
//             node_nbr[start_idx + j] = angle_nbr_pairs[j].second;
//         }
//     }
// }
// generate_vertex_neighbors(node_nbr_list2, numnbr2, cells, (int)cell_N / 3, N, nghst);
/**
 * @brief Generate vertex neighbors from triangular faces
 * @param faces Array of faces, each face contains 3 vertex indices
 * @param num_faces Number of faces
 * @param num_vertices Total number of vertices
 * @param num_neighbors Output array containing number of neighbors for each vertex
 * @param all_neighbors Output array containing all neighbors in flat format
 * @return Total number of neighbor relationships
 */
void generate_vertex_neighbors(int *all_neighbors, int *num_neighbors, int *faces, int num_faces, int num_vertices, int nghst){
    // int* faces, int num_faces, int num_vertices, int* num_neighbors, int* all_neighbors) {
    // Initialize vertex neighbor sets using vector<set<int>>
    vector<set<int>> vertex_neighbors(num_vertices);
    
    // For each face (triangle), add neighbor relationships
    for (int f = 0; f < num_faces; f++) {
        int v0 = faces[f * 3];
        int v1 = faces[f * 3 + 1]; 
        int v2 = faces[f * 3 + 2];
        
        // Each vertex in a triangle is a neighbor of the other two
        vertex_neighbors[v0].insert(v1);
        vertex_neighbors[v0].insert(v2);
        vertex_neighbors[v1].insert(v0);
        vertex_neighbors[v1].insert(v2);
        vertex_neighbors[v2].insert(v0);
        vertex_neighbors[v2].insert(v1);
    }
    
    // // Fill num_neighbors array and count total neighbors
    int total_neighbors = 0;
    for (int i = 0; i < num_vertices; i++) {
        num_neighbors[i] = static_cast<int>(vertex_neighbors[i].size());
        // total_neighbors += num_neighbors[i];
    }
    
    // Fill all_neighbors array with sorted neighbors for each vertex, pad with -1 if needed
    for (int i = 0; i < num_vertices; i++) {
        int j = 0;
        for (int neighbor : vertex_neighbors[i]) {
            all_neighbors[i * nghst + j] = neighbor;
            ++j;
        }
        // Pad the rest with -1
        for (; j < nghst; ++j) {
            all_neighbors[i * nghst + j] = -1;
        }
    }
}
/**
 * @brief Convert Python neighbours function to C++
 * This function processes simplices to find neighboring nodes and arranges them in anticlockwise order
 * Follows the exact Python logic: simpl is treated as having columns [r1, r2, r3] where r1 and r2 are the main edge vertices
 * @param Np Number of points
 * @param simpl Simplex array (3 columns: vertex1, vertex2, vertex3/-1)
 * @param vertices Array of vertex coordinates
 * @param nghst Maximum neighbors per vertex (default 12)
 * @return Pair of (cumlst, new_nbr) arrays
 */
pair<vector<int>, vector<int>> neighbours(int Np, vector<vector<int>> &simpl, Vec3d* vertices, 
    double boxL, int nghst){
    // Extract r1 and r2 from simpl (columns 0 and 1)
    vector<int> r1, r2;
    r1.reserve(simpl.size());
    r2.reserve(simpl.size());
    
    for (const auto &row : simpl) {
        if (row.size() >= 2) {
            r1.push_back(row[0]);
            r2.push_back(row[1]);
        }
    }

    int num_simpl = static_cast<int>(r1.size());
    
    // Calculate lst array - count neighbors for each point (Python: lst[i]=len(r1[r1==i])/2)
    vector<int> lst(Np, 0);
    for (int i = 0; i < Np; i++) {
        int count = 0;
        for (int j = 0; j < num_simpl; j++) {
            if (r1[j] == i) count++;
        }
        lst[i] = count / 2;  // Following Python logic exactly
    }
    
    // Calculate cumulative list (Python: cumlst[1:] = np.cumsum(lst))
    vector<int> cumlst(Np + 1, 0);
    for (int i = 1; i <= Np; i++) {
        cumlst[i] = cumlst[i-1] + lst[i-1];
    }
    
    // Create node_neighbour array (Python: node_neighbour[i]=r2[2*i])
    vector<int> node_neighbour(cumlst[Np]);
    int write_idx = 0;
    
    // Following the Python logic: for each vertex, collect its neighbors from r2 where r1 matches
    for (int i = 0; i < Np; i++) {
        for (int j = 0; j < num_simpl; j += 2) {  // Python uses 2*i indexing
            if (j < num_simpl && r1[j] == i) {
                if (write_idx < static_cast<int>(node_neighbour.size())) {
                    node_neighbour[write_idx] = r2[j];
                    write_idx++;
                }
            }
        }
    }
    
    // Sort neighbors in anticlockwise direction (new_way_nbrs equivalent)
    vector<int> new_nbr(Np * nghst, -1);
    
    for (int ip = 0; ip < Np; ip++) {
        // Get neighbors for current point
        vector<int> nbrs;
        for (int j = cumlst[ip]; j < cumlst[ip+1]; j++) {
            if (j < static_cast<int>(node_neighbour.size())) {
                nbrs.push_back(node_neighbour[j]);
            }
        }
        
        int num_nbr = static_cast<int>(nbrs.size());
        if (num_nbr == 0) continue;
        
        // Calculate angles for each neighbor
        vector<pair<double, int>> angle_nbr_pairs;
        for (int i = 0; i < num_nbr; i++) {
            int nbr_idx = nbrs[i];
            if (nbr_idx >= 0 && nbr_idx < Np) {
                double dx = vertices[nbr_idx].x - vertices[ip].x;
                if (dx >= 0.5 * boxL) dx -= boxL;
                if (dx < -0.5 * boxL) dx += boxL;
                double dy = vertices[nbr_idx].y - vertices[ip].y;
                if (dy >= 0.5 * boxL) dy -= boxL;
                if (dy < -0.5 * boxL) dy += boxL;
                double angle = atan2(dy, dx);
                if (angle < 0) {
                    angle += 2.0 * M_PI;
                }
                angle_nbr_pairs.emplace_back(angle, nbr_idx);
            }
        }

        // Sort by angle (anticlockwise)
        sort(angle_nbr_pairs.begin(), angle_nbr_pairs.end());
        
        // Store sorted neighbors
        int st_idx = ip * nghst;
        for (int j = 0; j < static_cast<int>(angle_nbr_pairs.size()) && (st_idx + j) < Np * nghst; j++) {
            new_nbr[st_idx + j] = angle_nbr_pairs[j].second;
        }
    }
    
    return make_pair(lst, new_nbr);
}
// /**
//  * @brief Overloaded version that works with C-style 3-column integer array
//  * @param Np Number of points
//  * @param simpl Simplex array as int* (3 columns: vertex1, vertex2, vertex3/-1)
//  * @param num_simpl Number of rows in simpl
//  * @param vertices Array of vertex coordinates
//  * @param nghst Maximum neighbors per vertex
//  * @return Pair of (cumlst, new_nbr) arrays
//  */
/**
 * @brief C-style interface that works with existing arrays
 * @param Np Number of points
 * @param simpl Simplex array (3 columns: vertex1, vertex2, vertex3/-1)
 * @param num_simpl Number of rows in simpl
 * @param vertices Array of vertex coordinates
 * @param cumlst Output cumulative list array (must be pre-allocated, size Np+1)
 * @param node_neighbour Output neighbor array (must be pre-allocated)
 * @param nghst Maximum neighbors per vertex
 */
void neighbours(int *numnbr, int *node_neighbour, vector<vector<int>> &simpl, Vec3d *vertices, 
    int Np, double boxL, int nghst = 12){
    auto result = neighbours(Np, simpl, vertices, boxL, nghst);
    // Copy cumlst
    for (size_t i = 0; i < result.first.size(); i++) {
        numnbr[i] = result.first[i];
    }
    
    // Copy new_nbr
    for (size_t i = 0; i < result.second.size(); i++) {
        node_neighbour[i] = result.second[i];
    }
}
//
bool is_wrapped_vertex(const Vec3d *pos, const vector<int> &nbrs, double boxL){
    double dx, dy, dz, dx_mi, dy_mi, dz_mi, d2_orig, d2_mi;
    for (size_t i = 0; i < nbrs.size(); i++){
        int n1 = nbrs[i];
        int n2 = nbrs[(i + 1) % nbrs.size()];
        dx = pos[n1].x - pos[n2].x;
        dy = pos[n1].y - pos[n2].y;
        dz = pos[n1].z - pos[n2].z;

        // Apply minimum image convention
        dx_mi = dx - boxL * round(dx / boxL);
        dy_mi = dy - boxL * round(dy / boxL);
        dz_mi = dz - boxL * round(dz / boxL);

        // If original != min-image, it’s a wrapped edge
        d2_orig = dx * dx + dy * dy + dz * dz;
        d2_mi = dx_mi * dx_mi + dy_mi * dy_mi + dz_mi * dz_mi;

        if (fabs(d2_orig - d2_mi) > 1e-8){
            return true;
        }
    }
    return false;
}
//
// using namespace MESH_p;
void determine_bdry_cdt(BoundaryType *btype, int *node_nbr_list, int *numnbr, 
        Vec3d *pos, const set<pair<int, int>> &edge_set, double boxL, int N, int nghst){
    // For each vertex, examine its neighbor ring. If an edge between consecutive neighbors
    // is missing from edge_set or is a wrapped edge (periodic wrap), mark as boundary.
    // If any boundary edge is wrapped -> PBC; else -> FBC. If no boundary edges -> NONE.
    for (int vi = 0; vi < N; ++vi){
        int num_nbr = numnbr[vi];
        if (num_nbr <= 0){
            btype[vi] = NONE;
            continue;
        }
        int start = vi * nghst;
        // collect neighbor indices
        int *nbradd = node_nbr_list + start;
        vector<int> nbrs(nbradd, nbradd+num_nbr);
        // Use the overload that inspects both topology (edge_set) and geometry (pos, boxL)
        bool is_bdry = is_boundary_vertex(nbrs, edge_set, pos, boxL);
        if (is_bdry){
            // If any of the boundary indicators are due to wrapping, mark PBC
            if (is_wrapped_vertex(pos, nbrs, boxL)){
                btype[vi] = PBC;
            } else {
                btype[vi] = FBC;
            }
        } else {
            btype[vi] = NONE;
        }
    }
}
#endif