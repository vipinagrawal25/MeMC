#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/Point_3.h>

#include <fstream>
#include <cassert>
#include <list>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <iostream>
#include <hdf5.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<K> GT;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<GT>       PDT;

typedef K::Point_3                                          Point_3;

typedef PDT::Vertex_handle                                  Vertex_handle;
typedef PDT::Point                                          Point;
typedef PDT::Iso_rectangle                                  Iso_rectangle;
typedef PDT::Finite_faces_iterator                          Face_iterator;

using namespace std;
// Function to sort simplices
vector<vector<int>> sort_simplices(const vector<vector<int>>& cells) {

   vector<vector<int>> nsimplices;
   for (const auto& scles : cells) {
      vector<int> nscles = scles;
      sort(nscles.begin(), nscles.end());
      nsimplices.push_back(nscles);

      // Generate permutations of the sorted simplex
      vector<int> temp(nscles);
      nsimplices.push_back({temp[1], temp[2], temp[0]});
      nsimplices.push_back({temp[2], temp[0], temp[1]});
      nsimplices.push_back({temp[0], temp[2], temp[1]});
      nsimplices.push_back({temp[1], temp[0], temp[2]});
      nsimplices.push_back({temp[2], temp[1], temp[0]});
   }

   // Sort the nsimplices
   sort(nsimplices.begin(), nsimplices.end(), [](const vector<int>& a, const vector<int>& b) {
      if (a[0] != b[0])
         return a[0] < b[0];
      else if(a[1]!=b[1])
         return a[1] < b[1];
      else
         return a[2] < b[2];
   });

   return nsimplices;
}

// Function to return cm_lst and node_nbr
pair<vector<int>, vector<int>> neighbours(int* lst, int Np,
   const vector<vector<int>>& simpl){
    
   vector<int> r1, r2, r3, cumlst(Np + 1, 0);
   for (const auto& s : simpl) {
      r1.push_back(s[0]);
      r2.push_back(s[1]);
      r3.push_back(s[2]);
   }

   for (int i = 0; i < Np; ++i) {
      lst[i] = count(r1.begin(), r1.end(), i) / 2;
   }

   partial_sum(lst, lst+Np, cumlst.begin() + 1);

   vector<int> node_neighbour(cumlst.back(), 0);
   for (size_t i = 0; i < node_neighbour.size(); ++i) {
      node_neighbour[i] = r2[2 * i];
   }

   return {cumlst, node_neighbour};
}

int hdf5_io_get_Np(string input_file){
   hid_t file_id,dataset_id;  /* identifiers */
   herr_t status;
   string dset_name="pos";

   if(access(input_file.c_str(),F_OK)!=0){
      fprintf(stderr, "The configuration file does not exit\n");
      exit(1);
   }

   file_id = H5Fopen(input_file.c_str(), H5F_ACC_RDONLY, 
      H5P_DEFAULT);
   dataset_id = H5Dopen(file_id, dset_name.c_str(),
      H5P_DEFAULT);

   hid_t dspace = H5Dget_space(dataset_id);
   const int ndims = H5Sget_simple_extent_ndims(dspace);
   hsize_t dims[ndims];
   H5Sget_simple_extent_dims(dspace, dims, NULL);

   return dims[0];
}

void hdf5_io_read_pos(double *Pos, string input_file){
   ///  @brief Read from the hdf5 file
   ///  @param Pos array containing co-ordinates of all the particles
   ///  @param input_file File name from which co-ordinate will be read

   hid_t file_id,dataset_id;  /* identifiers */
   herr_t status;
   string dset_name="pos";

   if(access(input_file.c_str(),F_OK)!=0){
      fprintf(stderr, "The configuration file does not exit\n");
      exit(1);
   }

   /* Open an existing file. */
   file_id = H5Fopen(input_file.c_str(), H5F_ACC_RDONLY, 
      H5P_DEFAULT);
   dataset_id = H5Dopen(file_id, dset_name.c_str(),
      H5P_DEFAULT);

   hid_t dspace = H5Dget_space(dataset_id);
   const int ndims = H5Sget_simple_extent_ndims(dspace);
   hsize_t dims[ndims];
   H5Sget_simple_extent_dims(dspace, dims, NULL);

   dataset_id = H5Dopen(file_id, dset_name.c_str(), 
      H5P_DEFAULT);
   status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, 
      H5S_ALL, H5P_DEFAULT, Pos);
   status = H5Dclose(dataset_id);
   status = H5Fclose(file_id);

   if(status != 0){
      fprintf(stderr, "file close failed\n");
   }
   
}

void cyclic_nbrs(int* new_nbr, const vector<int>& cmlist, 
                        const vector<int>& node_nbr,
                        Point* pts, int Np, int nghst, double lenth){

   // vector<int> new_nbr(nghst * Np, -1);

   double dx, dy;
   for (int ip = 0; ip < Np; ++ip){
      vector<int> nbrs(node_nbr.begin()+cmlist[ip], node_nbr.begin()+cmlist[ip + 1]);
      int num_nbr = (cmlist[ip + 1] - cmlist[ip]);
      vector<double> angles;

      for (int i : nbrs) {
         dx = pts[i].x() - pts[ip].x();
         if (dx >= 0.5 * lenth) dx -= lenth;
         if (dx < -0.5 * lenth) dx += lenth;

         dy = pts[i].y() - pts[ip].y();
         if (dy >= 0.5 * lenth) dy -= lenth;
         if (dy < -0.5 * lenth) dy += lenth;

         angles.push_back(atan2(dy, dx));
      }

      for (double& angle : angles){
         if (angle < 0)
            angle += 2 * M_PI;
      }

      vector<int> sorted(nbrs.size());
      iota(sorted.begin(), sorted.end(), 0);
      sort(sorted.begin(), sorted.end(), [&angles](int a, int b) {
         return angles[a] < angles[b];
      });

      int st_idx = ip * nghst;
      int end_idx = ip * nghst + num_nbr;

      for (int j = st_idx, k = 0; j < end_idx && k < nbrs.size(); ++j, ++k) {
         new_nbr[j] = nbrs[sorted[k]];
      }
   }
}

void hdf5_io_write_mesh(int *cmlist,
        int *node_nbr, int N, int ng, string output_file){

   ///  @brief Read the mesh from the hdf5 file
   ///  @param cmlist array containing the number of neighbours for each particle  
   ///  @param node_nbr array containing the list of neighbours for each particle  
   ///  @param input_file File name from which co-ordinate will be read

   hid_t   file_id, dset1, dataset_id, space_id;  /* identifiers */
   herr_t  status;
   int size_mesh; 
   hsize_t          dims; 

   size_mesh = ng*N;

   if(access(output_file.c_str(),F_OK)!=0){
      file_id = H5Fcreate (output_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   }else{
      file_id = H5Fopen (output_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
   }

   dims = N;
   space_id = H5Screate_simple (1, &dims, NULL);
   dset1 = H5Dcreate2(file_id, "/cumu_list", H5T_NATIVE_INT, space_id, H5P_DEFAULT,
            H5P_DEFAULT, H5P_DEFAULT);
   status = H5Dwrite (dset1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            cmlist);
   status = H5Dclose (dset1);
   status = H5Sclose (space_id);
   dims = size_mesh;
   space_id = H5Screate_simple (1, &dims, NULL);
   dset1 = H5Dcreate2(file_id, "/node_nbr", H5T_NATIVE_INT, space_id, H5P_DEFAULT,
            H5P_DEFAULT, H5P_DEFAULT);
   status = H5Dwrite (dset1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            node_nbr);
   status = H5Dclose (dset1);
   status = H5Sclose (space_id);
   status = H5Fclose(file_id);
   if(status != 0){
      fprintf(stderr, "file close failed\n");
   }
}

void hdf5_io_write_double(double *Pos, int N, 
        string input_file, string dset_name){

    ///  @brief hdf5 io for  the position 
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param N number of points 
    ///  @param input_file File name to dump all the co-ordinate
    /// 

    hid_t   file_id, dset1, space_id;  /* identifiers */
    herr_t  status;
    hsize_t          dims; 

  /* Open an existing file. */
    dims = N;
    file_id = H5Fcreate (input_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    space_id = H5Screate_simple (1, &dims, NULL);
    dset1 = H5Dcreate2(file_id, dset_name.c_str(), H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                Pos);
    status = H5Dclose (dset1);
    status = H5Sclose (space_id);
    status = H5Fclose (file_id);
  if(status != 0){
      fprintf(stderr, "file close failed\n");
  }
}

namespace memc{
   template <typename T>
   T* array(int nn, T value){
      T* arr;
      arr=(T*)calloc(nn,sizeof(T));
      for (int i = 0; i < nn; ++i){
         arr[i]=value;
      }
      return arr;
   }
}

int main(int argc, char const *argv[]){
   // vector<vector<int>> Pos;
   // Vec2d *Pos;
   Point* points;
   Point_3* pts_3d;
   int Np = (int)hdf5_io_get_Np(argv[1])/2;

   points = (Point *)calloc(Np, sizeof(Point));
   pts_3d = (Point_3*)calloc(Np, sizeof(Point_3));

   hdf5_io_read_pos((double *)points, argv[1]);

   double box_len=2*M_PI;
   int ng=12, k, cidx, k1;

   auto* nn_nbr=memc::array(ng*Np,-1);
   auto* ncml=memc::array(Np,0);

   Iso_rectangle domain(0, 0,box_len,box_len); // The cube for the periodic domain
   map<Vertex_handle, int> vertex_indices;
   vector<vector<int>>cells;

   PDT t(points,points+Np,domain);

   for(size_t i = 0; i < Np; ++i){
      Vertex_handle vh = t.nearest_vertex(points[i]);
      vertex_indices[vh] = i;
   }

   for (Face_iterator fit = t.finite_faces_begin(); fit != t.finite_faces_end(); ++fit){
      vector<int>tri;
      for (int i = 0; i < 3; ++i) { // Assuming 2D Delaunay triangulation (triangles)
         Vertex_handle vh = fit->vertex(i);
         tri.push_back(vertex_indices[vh]);
      }
      cells.push_back(tri);
   }

   vector<vector<int>> sorted_simplices = sort_simplices(cells);
   auto [cmlist, node_nbr] = neighbours(ncml, t.number_of_vertices(), sorted_simplices);
   cyclic_nbrs(nn_nbr, cmlist, node_nbr, points, Np, ng, box_len);

   for (int i = 0; i < Np; ++i){
      pts_3d[i] = Point_3(points[i].x(),points[i].y(),0.0);
      // cout << cmlist[i+1]-cmlist[i] - ncml[i] << "\t" << nn_nbr[i] << endl;
   }

   hdf5_io_write_double((double*)pts_3d, 3*Np, argv[2], "pos");
   hdf5_io_write_mesh(ncml, nn_nbr, Np, ng, argv[2]);

   free(ncml);
   free(nn_nbr);
   free(points);
   free(pts_3d);
   return 0;
}