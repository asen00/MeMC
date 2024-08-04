#include "metropolis.hpp"
#include "random_gen.hpp"
#include "multicomp.hpp"

#include <cmath>
#include <cstring>
#include <tuple>
#include <functional>
#include <algorithm>

const double pi = 3.14159265358979323846264;
// #include <cstdio>
// #include <iomanip>
// #include <sstream>
// #include <unistd.h>

// using namespace std;
extern "C" void  MC_listread(char *, double *, double *, bool *,int *, int *, 
                  bool *, int *, int *, double *, char *);
int get_nstart(int, int);
int McP::initMC(int N, string fname){
  char tmp_fname[128], temp_algo[128];
  string parafile, outfile;
  parafile = fname+"/para_file.in";
  sprintf(tmp_fname, "%s", parafile.c_str());
  MC_listread(temp_algo, &dfac, &kBT, &is_restart,
              &tot_mc_iter, &dump_skip, &is_fluid, &min_allowed_nbr,
              &fluidize_every, &fac_len_vertices, tmp_fname);
  algo=temp_algo;
  one_mc_iter = 2*N;
  dfac = sqrt(8*pi/(2*N-4))/dfac;
  acceptedmoves = 0;
  ofstream out_;
  out_.open( fname+"/mcpara.out");
  out_<< "# =========== monte carlo parameters ==========" << endl
      << " N " << N << endl
      << " algo = " << algo << endl
      << " dfac " << dfac << endl
      << " kbT " << kBT << endl
      << " is_restart " << is_restart << endl
      << " is_fluid " << is_fluid << endl
      << " tot_mc_iter " << tot_mc_iter << endl
      << " dump_skip " << dump_skip << endl
      << " min_allowed_nbr " << min_allowed_nbr << endl
      << " fluidize_every " << fluidize_every << endl;
  out_.close();
  return 1;
}
/*-----------------------*/
void McP::setEneVol() {
  EneMonitored = totEner;
  VolMonitored = totvol;
  // volt0 = totvol;
}
/*-----------------------*/
double McP::evalEnergy(Vec3d *Pos, MESH_p mesh){
  // if (fileptr.is_open()) {
  // fileptr << itr << "  " << (double)acceptedmoves/(double)one_mc_iter<< "  ";
  bende = beobj.bending_energy_total(Pos, mesh);
  stretche = steobj.stretch_energy_total(Pos, mesh);
  // It's okay to have if statement here.
  if (lipidobj.getcomp()>1) {regsole = lipidobj.reg_soln_tot(Pos, mesh);}
    // fileptr << bende << "  "<<stretche << "  ";
  totvol = steobj.volume_total(Pos, mesh);
    // totarea = steobj.area_total(Pos, mesh);
  // }
  totEner = bende+stretche+regsole;
  // cout << bende << " " << stretche << " " << regsole << endl;

  if (steobj.dopressure()){
    pre = -steobj.getpressure() * totvol;
    // fileptr << pre << "  ";
    totEner += pre;
  }
  
  // fileptr << totEner  << "  " << totvol  << endl;

  // if (is_fluid) {
  EneMonitored = totEner;
  VolMonitored = totvol;
  // }
  return totEner;
}
/*-----------------------*/
void McP::write_energy(fstream &fileptr, int itr){
  if (fileptr.is_open()){
    fileptr << itr << " " << (double)acceptedmoves/(double)one_mc_iter<< "  ";
    fileptr << bende << " " << stretche << "  ";
    if (lipidobj.getcomp()>1) fileptr << regsole << " ";
  }
  if (steobj.dopressure()) fileptr << pre << "  ";
  fileptr << EneMonitored  << "  " << VolMonitored  << endl;
}
/*-----------------------*/
// double McP::getarea(){return totarea;}
double McP::getvolume(){return totvol;}
bool McP::isrestart(){return is_restart;}
bool McP::isfluid(){return is_fluid;}
int McP::fluidizeevery(){return fluidize_every;}
int McP::dumpskip(){return dump_skip;}
int McP::totaliter(){return tot_mc_iter;}
int McP::onemciter(){return one_mc_iter;}
// int McP::monitoredVol(){return VolMonitored;}
// int McP::monitoredEn(){return EneMonitored;}
//
int del_nbr(int *nbrs, int numnbr, int idx){
  // delet int idx between i1 and i2 in the nbrs list
  int new_numnbr, delete_here;
  bool logic;

  logic = false;
  delete_here = 0;

  //  for (int i = 0; i < numnbr; ++i){
  //    if (nbrs[i]==nbrs[i+1]){
  //       // cout << "delnbr=" << idx << " "<< nbrs[delete_here] <<endl;
  //        for (int i = 0; i < numnbr; ++i){
  //           cout << nbrs[i] << " ";
  //        }
  //        cout << endl;
  //    }
  // }

  // for(int i=0; i<numnbr+3; i++)printf("%d \n", nbrs[i]);
  // printf("\n\n");

  while (!logic) {
    logic = (nbrs[delete_here] == idx);
    ++delete_here;
  }

  memcpy(nbrs + delete_here - 1, &nbrs[delete_here],
         sizeof(int) * (numnbr - delete_here + 1));

  // for(int i=0; i<numnbr+3; i++)printf("%d \n", nbrs[i]);
  // printf("\n\n");
  return numnbr - 1;
}

int add_nbr(int *nbrs, int numnbr, int idx, int i1, int i2) {
  // add int idx between i1 and i2 in the nbrs list
   int insert_here;
   bool logic;

   logic = false;
   insert_here = 0;

   // for(int i=0; i<numnbr+3; i++)printf("%d \n", nbrs[i]);
   //     printf("\n\n");

   // for (int i = 0; i < numnbr; ++i){
   //   if (nbrs[i]==nbrs[i+1]){
   //       // cout << "delnbr=" << idx << " "<< nbrs[delete_here] <<endl;
   //       for (int i = 0; i < numnbr; ++i){
   //          cout << nbrs[i] << " ";
   //       }
   //       cout << endl;
   //   }
   // }

   while (!logic) {
      logic = (nbrs[insert_here] == i1) || (nbrs[insert_here] == i2);
      ++insert_here;
   }

   logic = (nbrs[insert_here] == i1) || (nbrs[insert_here] == i2);
   if (logic) {
      memcpy(nbrs + insert_here, &nbrs[insert_here - 1],
           sizeof(int) * (numnbr - insert_here + 1));
      nbrs[insert_here] = idx;
   }else {
      insert_here = 0;
      memcpy(nbrs + insert_here, &nbrs[insert_here - 1],
           sizeof(int) * (numnbr - insert_here + 1));
      nbrs[insert_here] = idx;
   }

   // for(int i=0; i<numnbr+3; i++)printf("%d \n", nbrs[i]);
   //     printf("\n\n");

   return numnbr + 1;
}

bool McP::Boltzman(double DE, double activity){
   /// @brief Metropolis algorithm
   /// @param DE change in energy
   /// @param kbt boltzmann constant times temperature
   /// @return True if DE< 0 or the random number generated is less than
   /// exp(-DE/kbt)
   /// @details see
   /// https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm
   bool yes;
   double rand;
   DE += activity;
   yes = (DE <= 0.e0);
   if (!yes) {
      rand = RandomGenerator::generateUniform(0.0,1.0);
      yes = rand < exp(-DE / kBT);
   }
   return yes;
}
//
bool McP::Glauber(double DE, double activity){
  /// @brief Glauber algorithm
  /// @param DE change in energy
  /// @param kbt boltzmann constant times temperature
  bool yes;
  double rand;
  DE += activity;
  rand = RandomGenerator::generateUniform(0.0,1.0);
  yes = rand < 1 / (1 + exp(DE / kBT));
  return yes;
}
//
std::vector<double> McP::energy_mc_3d(Vec3d *pos, MESH_p mesh, int idx)
{
   double E_b, E_s, E_rs;
   vector<double> energy(3,0);
   int cm_idx, num_nbr;

   cm_idx = mesh.nghst * idx;
   num_nbr = mesh.numnbr[idx];

   E_b = beobj.bending_energy_ipart(pos, (int *)(mesh.node_nbr_list + cm_idx),
               num_nbr, idx, mesh.bdry_type, mesh.boxlen, mesh.edge);
   E_b += beobj.bending_energy_ipart_neighbour(pos, mesh, idx);
   E_s = steobj.stretch_energy_ipart(pos, (int *)(mesh.node_nbr_list + cm_idx),
               num_nbr, idx, mesh.nghst, mesh.bdry_type, mesh.boxlen, mesh.edge);

   energy[0]=E_b;
   energy[1]=E_s;
   // energy.push_back(E_b);
   // energy.push_back(E_s);
   // Try not to have if-else statements here.
   if (lipidobj.getcomp()>1){
      E_rs = lipidobj.gradphiener_ipart_all(pos, mesh, idx);
      energy[2] = E_rs;
   }
   //   if(st_p.do_stick)
   //   E_stick = lj_bottom_surface(pos[idx].z, st_p.is_attractive[idx],
   //       st_p.pos_bot_wall, st_p.epsilon, st_p.sigma); 

   //     if(afm.do_afm) E_afm = lj_afm(pos[idx], afm);

   //     if(spring.do_spring) E_spr = spring_energy(pos[idx], idx, mesh, spring);
   return energy;
}
//
int McP::monte_carlo_3d(Vec3d *pos, MESH_p mesh){
   int i, num_nbr, cm_idx;
   double x_o, y_o, z_o, x_n, y_n, z_n;
   double de, debe, dest, ders;
   // vector<double, double, double> Eini, Efin;
   double dxinc, dyinc, dzinc;
   double vol_i, vol_f;
   double dvol, de_vol, ini_vol, de_pressure;
   bool yes;
   int nframe;
  //
   nframe = get_nstart(mesh.N, mesh.bdry_type);
   acceptedmoves = 0;
   for (i = 0; i < one_mc_iter; i++) {
      int idx = RandomGenerator::intUniform(nframe, mesh.N-1);
      cm_idx = idx*mesh.nghst;
      num_nbr = mesh.numnbr[idx];
      auto Eini = energy_mc_3d(pos, mesh, idx);
      vol_i = steobj.volume_ipart(pos, (int *) (mesh.node_nbr_list + cm_idx),
               num_nbr, idx, mesh.bdry_type, mesh.boxlen, mesh.edge);
      //
      x_o = pos[idx].x; y_o = pos[idx].y; z_o = pos[idx].z;
      //
      dxinc = (dfac) * (RandomGenerator::generateUniform(-1.0,1.0));
      dyinc = (dfac) * (RandomGenerator::generateUniform(-1.0,1.0));
      dzinc = (dfac) * (RandomGenerator::generateUniform(-1.0,1.0));
      //
      x_n = x_o + dxinc; y_n = y_o + dyinc; z_n = z_o + dzinc;
       //
      pos[idx].x = x_n; pos[idx].y = y_n; pos[idx].z = z_n;
       //
      auto Efin = energy_mc_3d(pos, mesh, idx);

      debe=Efin[0]-Eini[0];
      dest=Efin[1]-Eini[1];
      ders=Efin[2]-Eini[2];

      de = debe+dest+ders;
      // cout << debe << "\t"  << dest << "\t" << ders << endl;
      //
      vol_f = steobj.volume_ipart(pos,
               (int *) (mesh.node_nbr_list + cm_idx), num_nbr, idx,
               mesh.bdry_type, mesh.boxlen, mesh.edge);
      dvol=0.5*(vol_f - vol_i);
      //
      yes = Boltzman(de, 0.0);
       // if(steobj.dovol()){
       // //   de_vol = vol_energy_change(mbrane, vol_p, dvol);
       // //   // cout << de << "\t";
       // //   de = de + de_vol;
       // //     
       // }

       // if(steobj.dopressure()){
       //   de_pressure = steobj.PV_change(dvol);
       //   de = de + de_pressure;
       // }

       // if (algo == "mpolis") {
       //   yes = Boltzman(de, 0.0);
       // } else if (algo == "glauber") {
       //   yes = Glauber(de, 0.0);
       // }
       //
      if(yes){
         acceptedmoves +=  1;
         EneMonitored += de;
         VolMonitored += 2*dvol;
         bende += debe;
         stretche += dest;
         regsole += ders;
         pre += de_pressure;
      } else {
         pos[idx].x = x_o;
         pos[idx].y = y_o;
         pos[idx].z = z_o;
      }
  }
  return acceptedmoves;
}
//
int McP::monte_carlo_fluid(Vec3d *pos, MESH_p mesh){

  int i, j, move;
  int nnbr_del1;
  int cm_idx_del1, cm_idx_del2;
  int cm_idx_add1, cm_idx_add2;
  int idx_del1, idx_del2;
  int idx_add1, idx_add2;
  int nframe;

  int nbr_add_1[12], nbr_add_2[12];
  int nbr_del_1[12], nbr_del_2[12];

  int N_nbr_del2, N_nbr_del1;
  int N_nbr_add2, N_nbr_add1;
  double det1, det2;
  Vec3d aft_ij;

  double KAPPA;
  double av_bond_len=mesh.av_bond_len;
  bool yes, logic;

  nframe = get_nstart(mesh.N, mesh.bdry_type);
  move = 0;

  int idxn, up, down;

  for (i = 0; i < one_mc_iter; i++) {
    // identify the pair to be divorced
    // stored as idx_del1 and idx_del2
    logic = false;
    while (!logic){
      idx_del1 = RandomGenerator::intUniform(nframe, mesh.N-1 );
      cm_idx_del1 = mesh.nghst * idx_del1;
      nnbr_del1 = mesh.numnbr[idx_del1];
      idxn = RandomGenerator::intUniform(0, mesh.nghst-1 );
      if (mesh.node_nbr_list[cm_idx_del1 + idxn] != -1) {
         idx_del2 = mesh.node_nbr_list[cm_idx_del1 + idxn];
         cm_idx_del2 = mesh.nghst * idx_del2;
         up = (idxn + 1 + nnbr_del1) % nnbr_del1;
         down = (idxn - 1 + nnbr_del1) % nnbr_del1;
         idx_add1 = mesh.node_nbr_list[idx_del1 * mesh.nghst + up];
         idx_add2 = mesh.node_nbr_list[idx_del1 * mesh.nghst + down];
         logic = idx_del2 > nframe && idx_add1 > nframe && idx_add2 > nframe;
      }else {
         logic = false;
      }
   }

   /* det1 = determinant(pos[idx_add1], pos[idx_add2], pos[idx_del1], mbrane.len); */
   /* det2 = determinant(pos[idx_add1], pos[idx_add2], pos[idx_del2], mbrane.len); */
   cm_idx_add1 = mesh.nghst * idx_add1;
   cm_idx_add2 = mesh.nghst * idx_add2;

   /* if (det1 * det2 < 0.0) { */
   aft_ij = pos[idx_add2] - pos[idx_add1];
   double dl = norm(aft_ij);
    // if (dl==0){
    //    cout << idx_add2 << " " << idx_add1 << endl;
    // }
    N_nbr_del1 = mesh.numnbr[idx_del1];
    N_nbr_del2 = mesh.numnbr[idx_del2];
    N_nbr_add1 = mesh.numnbr[idx_add1];
    N_nbr_add2 = mesh.numnbr[idx_add2];

   bool flip_condt1, flip_condt2, flip_condt3;
   bool accept_flip;

   flip_condt1 = (dl < fac_len_vertices*av_bond_len);
   flip_condt2 =  N_nbr_del1 > min_allowed_nbr && N_nbr_del2 > min_allowed_nbr;
   flip_condt3 =  N_nbr_add1 < 9 && N_nbr_add2 < 9;

   accept_flip = flip_condt1 && flip_condt2 && flip_condt3;

      if (accept_flip) {
         move = move + 1;
         /* print_sanity(pos, mesh.node_nbr_list+cm_idx_del1,
         * mesh.node_nbr_list+cm_idx_del2, */
         /*    mesh.node_nbr_list+cm_idx_add1,
         * mesh.node_nbr_list+cm_idx_add2, */
         /*         idx_del1, idx_del2, idx_add1, idx_add2, 
            (char*)"bef", i); */
         memcpy(nbr_del_1, &mesh.node_nbr_list[cm_idx_del1],
               sizeof(int) * mesh.nghst);
         memcpy(nbr_del_2, &mesh.node_nbr_list[cm_idx_del2],
               sizeof(int) * mesh.nghst);
         memcpy(nbr_add_1, &mesh.node_nbr_list[cm_idx_add1],
               sizeof(int) * mesh.nghst);
         memcpy(nbr_add_2, &mesh.node_nbr_list[cm_idx_add2],
               sizeof(int) * mesh.nghst);

        // form the bond
         N_nbr_add1 = add_nbr(nbr_add_1, mesh.numnbr[idx_add1], idx_add2,
                             idx_del1, idx_del2);
         N_nbr_add2 = add_nbr(nbr_add_2, mesh.numnbr[idx_add2], idx_add1,
                             idx_del1, idx_del2);

         // get divorced
         N_nbr_del1 = del_nbr(nbr_del_1, mesh.numnbr[idx_del1], idx_del2);
         N_nbr_del2 = del_nbr(nbr_del_2, mesh.numnbr[idx_del2], idx_del1);

         memcpy(mesh.node_nbr_list + cm_idx_del1, &nbr_del_1,
               sizeof(int) * mesh.nghst);
         memcpy(mesh.node_nbr_list + cm_idx_del2, &nbr_del_2,
               sizeof(int) * mesh.nghst);

         mesh.numnbr[idx_del1] = N_nbr_del1;
         mesh.numnbr[idx_del2] = N_nbr_del2;
         
         // for (int j= 0; j < mesh.numnbr[idx_del1]; ++j){
         //    if (mesh.node_nbr_list[cm_idx_del1+j] ==
         //        mesh.node_nbr_list[cm_idx_del1+j+1]){

         //       // cout << "idx_del1=" << idx_del1 << " " <<
         //       // nbr_del_1[j] << " " << nbr_del_1[j+1]
         //       // << endl;
         //    }
         // }

         // for (int j= 0; j < mesh.numnbr[idx_del2]; ++j){
         //    if (mesh.node_nbr_list[cm_idx_del2+j] ==
         //        mesh.node_nbr_list[cm_idx_del2+j+1]){
         //       // cout << "idx_del2="<< idx_del2 << " " <<
         //       // nbr_del_2[j] << " " << nbr_del_2[j+1]
         //       // << endl;
         //    }
         // }

         memcpy(mesh.node_nbr_list + cm_idx_add1, &nbr_add_1,
               sizeof(int) * mesh.nghst);
         memcpy(mesh.node_nbr_list + cm_idx_add2, &nbr_add_2,
               sizeof(int) * mesh.nghst);

         mesh.numnbr[idx_add1] = N_nbr_add1;
         mesh.numnbr[idx_add2] = N_nbr_add2;

         // for (int j= 0; j < mesh.numnbr[idx_add1]; ++j){
         //    if (mesh.node_nbr_list[cm_idx_add1+j] ==
         //        mesh.node_nbr_list[cm_idx_add1+j+1]){
         //       // cout<< "idx_add1=" << idx_add1 << " "
         //       // << nbr_add_1[j] << " " << nbr_add_1[j+1]
         //       // << endl;
         //    }
         // }

         // for (int j= 0; j < mesh.numnbr[idx_add2]; ++j){
         //    if (mesh.node_nbr_list[cm_idx_add2+j] ==
         //        mesh.node_nbr_list[cm_idx_add2+j+1]){
         //       // cout<< "idx_add2=" << idx_add2 << " "
         //       // << nbr_add_2[j] << " " << nbr_add_2[j+1]
         //       // << endl;
         //    }
         // }


      }
    /* } */
  }
  return move;
}
//
int McP::monte_carlo_lipid(Vec3d *pos, MESH_p mesh){
   int exchngdmoves = 0;
   int idx1, idx2, cm_idx1;
   int nframe = get_nstart(mesh.N, mesh.bdry_type);
   double Eini, Efin;
   bool yes, logic;
   int lip_idx1, lip_idx2, idxn, logic_break;
   for (int i = 0; i < one_mc_iter; ++i){
      logic = true;
      idx1 = RandomGenerator::intUniform(nframe, mesh.N-1);
      cm_idx1 = mesh.nghst * idx1;

      idxn=0;
      logic_break=0;
      while(logic && logic_break<2*mesh.nghst){
         idxn = RandomGenerator::intUniform(0, mesh.nghst-1);
         idx2 = mesh.node_nbr_list[cm_idx1+idxn];
         if (idx2!=-1){
            logic = lipidobj.lipA[idx1] == lipidobj.lipA[idx2];
         }
         ++logic_break;
      }

      // idxn=RandomGenerator::intUniform(0, mesh.nghst-1);
      // idx2 = mesh.node_nbr_list[cm_idx1+idxn];
      // if (idx2!=-1){
      //    logic = lipidobj.lipA[idx1] == lipidobj.lipA[idx2];
      // }

      if (!logic){
         lip_idx1 = lipidobj.lipA[idx1];
         lip_idx2 = lipidobj.lipA[idx2];

         Eini = lipidobj.reg_soln_ipart(pos, mesh, idx1).x;
         Eini += lipidobj.reg_soln_ipart_neighbour(pos, mesh, idx1);
         Eini += lipidobj.reg_soln_ipart(pos, mesh, idx2).x;
         Eini += lipidobj.reg_soln_ipart_neighbour(pos, mesh, idx2);

         lipidobj.lipA[idx2] = lip_idx1;
         lipidobj.lipA[idx1] = lip_idx2;

         Efin = lipidobj.reg_soln_ipart(pos, mesh, idx1).x;
         Efin += lipidobj.reg_soln_ipart_neighbour(pos, mesh, idx1);
         Efin += lipidobj.reg_soln_ipart(pos, mesh, idx2).x;
         Efin += lipidobj.reg_soln_ipart_neighbour(pos, mesh, idx2);

         yes = Boltzman(Efin-Eini, 0.0);
         // cout << exp(-Efin+Eini) << endl;
         if (yes){
            ++exchngdmoves;
         }else{
            lipidobj.lipA[idx1] = lip_idx1;
            lipidobj.lipA[idx2] = lip_idx2;
         }
      }
   }
   return exchngdmoves;
}