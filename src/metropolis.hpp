#ifndef METROPOLIS_HPP
#define METROPOLIS_HPP
#include <string>
#include <fstream>
#include <vector>

#include "mesh.hpp"
#include "vector.hpp"
#include "bending.hpp"
#include "stretching.hpp"
#include "multicomp.hpp"
class McP{
public : 
  McP (BE &beobj, STE &steobj, MulCom &lipidobj): 
  beobj(beobj), steobj(steobj), lipidobj(lipidobj) {};
  int monte_carlo_3d(Vec3d *pos, MESH_p mesh);
  std::vector<double> energy_mc_3d(Vec3d *pos, MESH_p mesh, int);
  int monte_carlo_fluid(Vec3d *, MESH_p);
  int monte_carlo_lipid(Vec3d *pos, MESH_p mesh);
  bool Boltzman(double DE, double activity);
  bool Glauber(double DE, double activity);
  int initMC(int, std::string);
  bool isfluid();
  bool isrestart();
  int fluidizeevery();
  int dumpskip();
  int totaliter();
  int onemciter();
  double evalEnergy(Vec3d *, MESH_p);
  void write_energy(fstream &, int);
  double getarea();
  double getvolume();
  void setEneVol();
private:
  BE &beobj;
  STE &steobj;
  MulCom &lipidobj;
  std::string algo;
  double dfac;
  int one_mc_iter, tot_mc_iter, dump_skip;
  double kBT;
  double delta; // increment of position
  bool is_restart;
  bool is_fluid;
  int min_allowed_nbr;
  int fluidize_every;
  double fac_len_vertices;
  double totEner, totvol, bende, stretche, pre, regsole;
  double EneMonitored, VolMonitored;
  double volt0;
  int acceptedmoves;
};
#endif