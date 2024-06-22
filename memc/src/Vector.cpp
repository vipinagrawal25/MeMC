#include "Vector.hpp"

double inner_product(Vec3d s1, Vec3d s2){
    return s1.x*s2.x + s1.y*s2.y + s1.z*s2.z;
}

 double norm(Vec3d s1){
    return sqrt(inner_product(s1,s1));
}

 double normsq(Vec3d s1){
    return inner_product(s1,s1);
}

Vec3d Vec3d_add(Vec3d s1, Vec3d s2, double fac){
    /* Returns s1 + fac*s2*/
    Vec3d add;
    add.x = s1.x + fac*s2.x;
    add.y = s1.y + fac*s2.y;
    add.z = s1.z + fac*s2.z;
    return add;
}

Vec3d cross_product(Vec3d s1, Vec3d s2){

    Vec3d crosprod;
    crosprod.x = s1.y*s2.z - s1.z*s2.y;
    crosprod.y = s1.z*s2.x - s1.x*s2.z;
    crosprod.z = s1.x*s2.y - s1.y*s2.x;
    return crosprod;
}
/*---------------------------------------*/
Vec3d diff_pbc(Vec3d r1, Vec3d r2, double len){
    Vec3d rij;

    rij = r2 - r1;

#ifdef flat
    if(rij.x >= 0.5*len) rij.x = rij.x - len;
    if(rij.y >= 0.5*len) rij.y = rij.y - len;

    if(rij.x < -0.5*len) rij.x = rij.x + len;
    if(rij.y < -0.5*len) rij.y = rij.y + len;
#endif

    return rij;
}

