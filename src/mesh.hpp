#ifndef MESH_HPP
#define MESH_HPP
#include <string>
typedef struct{
    /// @brief Mesh Structure
    /// @param numnbr; number of neighbours
    /// @param node_nbr_list; list of neighbours of a node
    int N, bdry_type;
    int nghst;
    std::string topology;
    int *numnbr;
    int *node_nbr_list;
    double boxlen;
    int edge; // storing corner index specially for periodic case.
    double av_bond_len;
}MESH_p;
extern "C" void MeshRead(int *, int *, double *, char *);
#endif