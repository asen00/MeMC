#ifndef MULTICOMP_HPP
#define MULTICOMP_HPP
#include "misc.hpp"
#include "mesh.hpp"
#include "vector.hpp"
class MulCom{
public:
  Vec2d reg_soln_ipart(Vec3d *pos, MESH_p mesh, int idx);
  double reg_soln_ipart_neighbour(Vec3d *pos, MESH_p mesh, int idx);
  double reg_soln_tot(Vec3d *pos, MESH_p mesh);
private:
  bool islipid;
  bool *lipA;
  double phi_ipart(int *node_nbr, int num_nbr, int idx);
  void phi_ipart_neighbour(double *phi, MESH_p mesh, int idx);
  Vec2d gradphi_sq(double *phi, Vec3d *pos, int *node_nbr, int num_nbr, int idx,
    int bdry_type, double lenth, int edge);
};
#endif