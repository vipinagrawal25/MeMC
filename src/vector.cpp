#include "vector.hpp"

Vec3d Vec3d_add(Vec3d s1, Vec3d s2, double fac){
    /* Returns s1 + fac*s2*/
    Vec3d add;
    add.x = s1.x + fac*s2.x;
    add.y = s1.y + fac*s2.y;
    add.z = s1.z + fac*s2.z;
    return add;
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