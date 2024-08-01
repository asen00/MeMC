#include <ctime>
#include <algorithm>
#include <fstream>

#include "multicomp.hpp"
#include "vector.hpp"

extern "C" void LipidRead(int *, double *, double *, char *);

void fillPoints(std::vector<int>& points, double fraction, int N){
   
   int numOnes = static_cast<int>(fraction * N);
   int numZeros = N - numOnes;

   // Fill the vector with the required number of 1s and 0s
   for (int i = 0; i < numOnes; ++i) points.push_back(1);
   for (int i = numOnes; i < N; ++i) points.push_back(0);

   // Shuffle the vector to mix 0s and 1s
   std::srand(unsigned(std::time(0)));
   std::random_shuffle(points.begin(), points.end());
}
//
int MulCom::initMulCom(int N, double av_bond_len, BE &bendobj, std::string fname){
    char tmp_fname[128];
    string parafile, outfile;

    parafile = fname+"/para_file.in";
    sprintf(tmp_fname, "%s", parafile.c_str() );

    LipidRead(&ncomp, &kai, &lipfrac, tmp_fname);
    epssqby2=(4-2*kai)*av_bond_len*av_bond_len;

    fillPoints(lipA, lipfrac, N);
   
    bendobj.init_coefbend(lipA, N);

    std::ofstream out_;
    out_.open(fname+"/lipidpara.out");
    out_<< "# =========== Lipid parameters ==========" << endl
      << " N " << N << endl
      << " ncom = " << ncomp << endl
      << " kai " << kai << endl
      << " lipfrac " << lipfrac << endl
      << " epssqby2 = " << epssqby2 << endl;
    out_.close();

    out_.open(fname+"/lipA.txt");
    for (int value : lipA){
        out_ << value << "\n";  // Save each element on a new line
    }
    out_.close();

    return 0;
}
//
Vec2d MulCom::gradphi_sq(double *phi, Vec3d *pos, int *node_nbr, int num_nbr, 
    int idx, int bdry_type, double lenth, int edge){
    /// @brief Per-vertex gradient estimation in DOI: 10.2312/stag.20181301
    int j, jdx, kdx, k;
    Vec3d  rk, ri, rj, rkp;
    Vec3d rji, ip1, rik;
    double area1 = 0e0;
    ri = pos[idx];
    Vec3d gradphi;
    Vec2d grad_ar;
    if (bdry_type==1||idx > edge){
        for (j = 0; j < num_nbr; j++){
            jdx = node_nbr[j];
            k = (j+1)%num_nbr;
            kdx = node_nbr[k];
            rj = pos[jdx]; rk = pos[kdx];
            rji  = rj - ri;
            rik  = ri - rk;
            ip1 = cross_product(rik, rji);
            area1 = area1 + 0.5*norm(ip1);
            ip1 = ip1*1/norm(ip1);
            gradphi=gradphi+(cross_product(rik, ip1)/2)*(phi[j+1]-phi[0])
                    + (cross_product(rji, ip1)/2)*(phi[k+1]-phi[0]);
        }
    }else{
        for (j = 0; j < num_nbr; j++){
            jdx = node_nbr[j];
            k = (j+1)%num_nbr;
            kdx = node_nbr[k];
            rj = pos[jdx]; rk = pos[kdx];
            rji = change_vec_pbc(rj-ri,lenth);
            rik = change_vec_pbc(ri-rk,lenth);
            ip1 = cross_product(rik, rji);
            area1 = area1 + 0.5*norm(ip1);
            ip1 = ip1*1/norm(ip1);
            gradphi=gradphi+(cross_product(rik, ip1)/2)*(phi[j+1]-phi[0])
                    + (cross_product(rji, ip1)/2)*(phi[k+1]-phi[0]);
        }
    }
    gradphi = gradphi/area1;
    grad_ar.x = normsq(gradphi);
    grad_ar.y = area1;
    return grad_ar;
}
//
double MulCom::phi_ipart(int *node_nbr, int num_nbr, int idx){
    double phi = (double) lipA[idx];
    int jdx;
    for (int i =0; i < num_nbr; i++){
        jdx = node_nbr[i];
        phi+= (double) lipA[jdx];
    }
    return phi/(1+num_nbr);
}
//
void MulCom::phi_ipart_neighbour(double *phi, MESH_p mesh, int idx){
    int cnt,j;
    int num_nbr_j;
    int nbr, cm_idx_nbr;
    cnt = 0;
    int cm_idx = mesh.nghst * idx;
    int num_nbr = mesh.numnbr[idx];
    phi[cnt] = phi_ipart((int *)(mesh.node_nbr_list + cm_idx), num_nbr, idx);
    for (j = idx*mesh.nghst; j < idx*mesh.nghst + mesh.numnbr[idx]; j++){
        cnt+=1;
        nbr = mesh.node_nbr_list[j];
        num_nbr_j = mesh.numnbr[nbr];
        cm_idx_nbr = nbr*mesh.nghst;
        phi[cnt] = phi_ipart((int *) mesh.node_nbr_list + cm_idx_nbr, num_nbr_j,
                    nbr);
    }
}
//
Vec2d MulCom::reg_soln_ipart(Vec3d *pos, MESH_p mesh, int idx){
    /// TODO: Add interaction (kai*phi*(1-phi))
    int num_nbr = mesh.numnbr[idx];
    double phi[1+num_nbr];
    phi_ipart_neighbour(phi, mesh, idx);
    Vec2d grad_arr, out_array(0.0,0.0);
    if(phi[0]==0||phi[0]==1) return out_array;
    double mixenergy = phi[0]*log(phi[0]) + (1-phi[0])*log(1-phi[0])
                        + kai*phi[0]*(1-phi[0]);
    // phi has a size of 1+num_nbr. You can't access phi[idx]
    int cm_idx = idx*mesh.nghst;
    grad_arr = gradphi_sq(phi, pos, (int *)(mesh.node_nbr_list + cm_idx),
                num_nbr, idx, mesh.bdry_type, mesh.boxlen, mesh.edge);
    out_array.x = mixenergy + epssqby2*grad_arr.x; /// actual free energy
    out_array.y = grad_arr.y;
    return out_array;
}
//
double MulCom::reg_soln_ipart_neighbour(Vec3d *pos, MESH_p mesh, int idx){
    int j, nbr;
    double regsol_fe = 0e0;
    Vec2d regsol_ar;
    for (j = idx*mesh.nghst; j < idx*mesh.nghst + mesh.numnbr[idx]; j++){
        nbr = mesh.node_nbr_list[j];
        regsol_ar = reg_soln_ipart(pos, mesh, nbr);
        regsol_fe +=  regsol_ar.x;
   }
   return regsol_fe;
}
//
double MulCom::reg_soln_tot(Vec3d *pos, MESH_p mesh){
    int idx, st_idx;
    int num_nbr, cm_idx;
    double rsfe_tot;
    Vec2d rs_ar;
    rsfe_tot = 0e0;
    st_idx = get_nstart(mesh.N, mesh.bdry_type);
    for(idx = st_idx; idx < mesh.N; idx++){
        cm_idx = idx*mesh.nghst;
        num_nbr = mesh.numnbr[idx];
        rs_ar=reg_soln_ipart(pos, mesh, idx);
        rsfe_tot = rsfe_tot + rs_ar.x;
    }
    return rsfe_tot;
}