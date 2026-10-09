#ifndef MESH_HPP
#define MESH_HPP
typedef struct{
    /// @brief Mesh Structure
    /// @param numnbr; number of neighbours
    /// @param node_nbr_list; list of neighbours of a node
    int N,  bdry_type;
    int nghst;
    double len_x;
    double len_y; 
    int *numnbr;
    int *node_nbr_list;
}MESH_p;
//
#endif
