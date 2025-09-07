#include "mesh.hpp"
#include "vector.hpp"
#include <vector>
#include <set>
#include <utility>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <random>
#include <preprocess/flat.hpp>
#include "misc.hpp"

using namespace std;
bool is_boundary_vertex(const vector<int> &nbrs, const set<pair<int, int>> &edge_set);
void fillPoints(vector<int> &points, double fraction, int N);
void fillPoints(int *points, double fraction, int N);
void fillPoints(int *points, int N, int value);
vector<vector<int>> sort_simplices(int *cells, int num_cells);
void get_neighbours(const vector<vector<int>> &simplices, int *numnbr,
                    int *node_nbr, int N, int nghst);
vector<vector<int>> remove_duplicates(vector<vector<int>> &sorted_cells);
// void MESH::initMESH(){
//     MeshRead(&bdry_type, &nghst, &radius, tmp_fname);
//     N = (int)hdf5_io_get_Np(outfolder+"/input.h5", "pos")/3;
//     Pos = (Vec3d *)calloc(N, sizeof(Vec3d));
//     numnbr = (int *)calloc(N, sizeof(int));
//     node_nbr_list = (int *)calloc(nghst*N, sizeof(int));
// }

pair<double, double> MESH_p::get_box_dim()
{
    // Vec3d *Pos = mesh.pos;
    double xmin = 1e7, xmax = -1e7, ymin = 1e7, ymax = -1e7;
    for (int i = 0; i < N; ++i){
        if (pos[i].x < xmin)
            xmin = pos[i].x;
        if (pos[i].x > xmax)
            xmax = pos[i].x;
        if (pos[i].y < ymin)
            ymin = pos[i].y;
        if (pos[i].y > ymax)
            ymax = pos[i].y;
    }
    double xlen = xmax - xmin;
    double ylen = ymax - ymin;
    return {xlen, ylen};
}

bool MESH_p::isPlaner() {
    for (int i = 0; i < N; i++)
    {
        if (pos[i].z != 0)
        {
            return false;
        }
    }
    return true;
}

double MESH_p::calculateRadius(){
    double sum_distances = 0.0;

    for (int i = 0; i < N; i++)
    {
        auto point = pos[i];
        double distance = sqrt(point.x * point.x + point.y * point.y + point.z * point.z);
        sum_distances += distance;
    }

    // Calculate average distance (approximate radius)
    double radius = sum_distances / N;
    return radius;
}

bool MESH_p::determine_pbc(){
    // Check if any node has a neighbor that wraps around the boundary
    for (int i = 0; i < N; ++i)
    {
        int current_index = i * nghst;
        for (int j = 0; j < numnbr[i]; ++j)
        {
            int neighbor = node_nbr_list[current_index + j];
            if (abs(pos[i].x - pos[neighbor].x) > boxlen / 2 ||
                abs(pos[i].y - pos[neighbor].y) > boxlen / 2 ||
                abs(pos[i].z - pos[neighbor].z) > boxlen / 2)
            {
                return true;
            }
        }
    }
    return false;
}

MESH_p::MESH_p(string outfolder){
    char tmp_fname[128], tmp_dist[128];
    string para_file = outfolder + "/para_file.in";
    sprintf(tmp_fname, "%s", para_file.c_str());

    MeshRead(&bdry_type, &compfrac, tmp_dist, tmp_fname);
    distribution = tmp_dist;
    N = (int) hdf5_io_get_Np(outfolder + "/input.h5", "pos") / 3;
    if (compfrac > 1 || compfrac < 0){
        cerr << "Error: The fraction of the component should be between 0 and 1" << endl;
        exit(EXIT_FAILURE);
    }

    if (distribution != "Random" && distribution != "random" &&
        distribution != "Janus" && distribution != "janus" &&
        distribution != "Point" && distribution != "point"){
        cerr << "Error: The distribution type is not recognized. Please use either Random, Janus or Point." << endl;
        exit(EXIT_FAILURE);
    }

    if (compfrac > 0 && compfrac <= 1){
        ncomp = 2;
    }
    else{
        ncomp = 1;
    }

    if (bdry_type != 0 && bdry_type != 1 && bdry_type != 2){
        cerr << "Error: The boundary type is not recognized. Please use either 0, 1 or 2." << endl;
        exit(EXIT_FAILURE);
    }

    pos = new Vec3d[N];
    numnbr = new int[N];
    node_nbr_list = new int[N * nghst];

    int *numnbr2 = new int[N];
    int *node_nbr_list2 = new int[N * nghst];

    compA = new int[N];

    hdf5_io_read_double((double *)pos, outfolder + "/input.h5", "pos");
    if (isPlaner()){
        sphere = false;
        //
        if (hdf5_io_has_dataset(outfolder + "/input.h5", "cumu_list") &&
            hdf5_io_has_dataset(outfolder + "/input.h5", "node_nbr")){
            // If both "cumu_list" and "node_nbr" exist, read mesh connectivity data.
            hdf5_io_read_mesh((int *)numnbr, (int *)node_nbr_list, outfolder + "/input.h5");
        }else if (hdf5_io_has_dataset(outfolder + "/input.h5", "cells")){
            auto cell_N = (int)hdf5_io_get_Np(outfolder + "/input.h5", "cells");
            cout << "Number of cells (triangles) in the mesh: " << (int)cell_N / 3 << endl;
            cells = new int[cell_N]; // Allocate memory for cells (6N points, each with x,y,z)
            hdf5_io_read_int((int *)cells, outfolder + "/input.h5", "cells");
            auto sort_tri = sort_simplices(cells, (int)cell_N / 3);
            //
            // The next step is important because you do not know if you are getting faces directly after delaunay triangulation or
            // preprocessed such that triangles are stored for every vertex.
            vector<vector<int>> unique_tri = remove_duplicates(sort_tri);
            get_neighbours(node_nbr_list,  numnbr, unique_tri, N, nghst);
        }else{
            cerr << "Error: The input HDF5 file must contain either 'cumu_list' and 'node_nbr' datasets or a 'cells' dataset." << endl;
            cerr << "First compute the Delaunay triangulation of the system and rerun the code." << endl;
            exit(EXIT_FAILURE);
        }
        //
        boxlen = get_box_dim().first * (1 + 1 / sqrt(N));
        sort_nbrs(node_nbr_list, pos, N, numnbr, nghst, boxlen);
        lastbdry = put_boundary_first(pos, node_nbr_list, numnbr, N, nghst, boxlen);
        cout<< "Last boundary index after rearranging: " << lastbdry << endl;
        //
        if (lastbdry == N - 1){
            cout << "Warning: All the points are boundary points. The mesh is likely to be incorrect." << endl;
        }
        if (lastbdry == -1){
            cout << "Warning: No boundary points found. The mesh is likely to be incorrect." << endl;
        }
        // Determine if the mesh is periodic or not.
        pbc = determine_pbc();
        if (pbc){
            if (bdry_type == 0 || bdry_type == 1){
                cerr << "Error: Boundary type of fixed frame (0) or channel (1) cannot be used with periodic mesh.\n"
                        "Run this code with the right boundary condition"
                     << endl;
                exit(EXIT_FAILURE);
            }
        }else{
            order_boundary_neighbors(node_nbr_list, numnbr, pos, lastbdry, nghst);
        }
    }
    else{
        sphere = true;
        lastbdry = -1;
        boxlen = 0;
        bdry_type = 2; // Sphere always has a free boundary condition.
        radius = calculateRadius();
        cout << "radius = " << radius << endl;
        ini_vol = 4e0 / 3e0 * M_PI * radius * radius * radius;
        hdf5_io_read_mesh((int *)numnbr, (int *)node_nbr_list, outfolder + "/input.h5");
    }

    if (distribution == "Random" || distribution == "random"){
        fillPoints(compA, compfrac, N);
    }
    else if (distribution == "Janus" || distribution == "janus"){
        zattr = 2 * (compfrac - 0.5) * radius;
        identify_attractive_part(compA, pos, zattr, N);
    }
    else if (distribution == "Point" || distribution == "point"){
        compA[0] = 1;
    }

    for (i = 0; i < N; i++){
        num_nbr = numnbr[i];
        cm_idx = nghst * i;
        for (k = cm_idx; k < cm_idx + num_nbr; k++)
        {
            j = node_nbr_list[k];
            dr = diff_pbc(pos[j], pos[i], boxlen);
            sum_lij += sqrt(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);
            npairs++;
        }
    }
    av_bond_len = sum_lij / npairs;
}

void MESH_p::free(){
    // @brief This is something similar to destructor of the class.
    // I do not remember now but there was some issue in defining the normal
    // destructor of the structure
    delete[] pos;
    delete[] cells;  // Clean up cells memory
    delete[] node_nbr_list;
    delete[] numnbr;
    delete[] compA;
}

void fillPoints(vector<int> &points, double fraction, int N){
    int numOnes = static_cast<int>(fraction * N);
    int numZeros = N - numOnes;

    // Fill the vector with the required number of 1s and 0s
    points.clear(); // Clear the vector first if you're reusing it
    for (int i = 0; i < numOnes; ++i)
        points.push_back(1);
    for (int i = numOnes; i < N; ++i)
        points.push_back(0);

    // Shuffle the vector using shuffle
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    shuffle(points.begin(), points.end(), default_random_engine(seed));
}

void fillPoints(int *points, double fraction, int N){
    int numOnes = static_cast<int>(fraction * N);
    int numZeros = N - numOnes;

    // Fill the array with the required number of 1s and 0s
    for (int i = 0; i < numOnes; ++i)
        points[i] = 1;
    for (int i = numOnes; i < N; ++i)
        points[i] = 0;

    // Shuffle the array
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(points, points + N, std::default_random_engine(seed));
}

void fillPoints(int *points, int N, int value){
    // Fill the array with the required number of 1s and 0s
    for (int i = 0; i < N; ++i)
        points[i] = value;
}

// int MESH_p::get_bdry(const set<pair<int, int>> &edge_set){
//     /// @brief This function assumes that the boundary vertices are at the beginning of the list.
//     /// @brief This function returns the same boundary points irrespective of the periodicity or
//     ///                aperiodicity of the mesh.
//     /// @param edge_set Set of edges in the mesh
//     for (int i = 0; i < N; ++i)
//     {
//         int *nbrs = node_nbr_list + i * nghst;
//         vector<int> nbr_list(nbrs, nbrs + numnbr[i]);
//         if (is_boundary_vertex(nbr_list, edge_set) == false)
//             return i - 1;
//     }
//     return -1;
// }

// void get_neighbours(const vector<vector<int>> &simplices, int *numnbr, int *node_nbr, int N){
//     vector<int> cumlst(N+1, 0);
//     // First count occurrences of each vertex as first vertex
//     vector<int> vertex_count(N, 0);
//     for (const auto &simplex : simplices) {
//         vertex_count[simplex[0]]++;
//     }

//     // Calculate cumulative sum for cumlst
//     cumlst[0] = 0;
//     for (int i = 0; i < N; i++) {
//         cumlst[i + 1] = cumlst[i] + vertex_count[i]/2;  // Divide by 2 as in Python
//     }

//     // Fill node_neighbour array
//     vector<int> current_pos(N, 0);  // Track current position for each vertex
//     for (size_t i = 0; i < simplices.size(); i += 2) {  // Step by 2 as we only want half the duplicates
//         int v1 = simplices[i][0];    // First vertex
//         int v2 = simplices[i][1];    // Second vertex (neighbor)
        
//         // Add v2 as neighbor of v1
//         int pos = cumlst[v1] + current_pos[v1];
//         node_nbr[pos] = v2;
//         current_pos[v1]++;
//     }

//     for (int i = 0; i < N; ++i) {
//         numnbr[i] = cumlst[i + 1] - cumlst[i];
//     }
// }

