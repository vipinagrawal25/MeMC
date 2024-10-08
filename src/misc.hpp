#ifndef MISC_HPP
#define MISC_HPP
int get_nstart(int N, int bdrytype);
#include "vector.hpp"
#include <vector>

void print(int *arr, int nn);
void print(double *arr, int nn);
void fillPoints(std::vector<int>& points, double fraction, int N);
void fillPoints(int *points, double fraction, int N);
void fillPoints(int *points, int N, int value);

template<typename T>
void print(const std::vector<T>& vec){
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i] << " ";
    }
    std::cout << std::endl;
}
/*-----------------------------------------------*/
inline Vec3d change_vec_pbc(Vec3d ab, double lenth){
   if (ab.x >= 0.5 * lenth) ab.x -= lenth;
   if (ab.x < -0.5 * lenth) ab.x += lenth;

   if (ab.y >= 0.5 * lenth) ab.y -= lenth;
   if (ab.y < -0.5 * lenth) ab.y += lenth;

   return ab;
}

inline Vec3d diff_pbc(Vec3d a, Vec3d b, double lenth){
   return change_vec_pbc(a-b, lenth);
}

#endif