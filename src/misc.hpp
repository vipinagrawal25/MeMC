#ifndef MISC_HPP
#define MISC_HPP
int get_nstart(int N, int bdrytype);
#include "vector.hpp"
#include <vector>

void print(int *arr, int nn);
void print(double *arr, int nn);
void identify_attractive_part(int *is_attractive, Vec3d *pos,
    double theta_attr, int N);

template<typename T>
void print(const std::vector<T>& vec){
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i] << " ";
    }
    std::cout << std::endl;
}
// /*-----------------------------------------------*/

#endif