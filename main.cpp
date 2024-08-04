#include <iostream>
#include <cstdlib>
#include <time.h>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <stdbool.h>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <fstream>
#include <iostream>

#include "metropolis.hpp"
#include "bending.hpp"
#include "stretching.hpp"
#include "random_gen.hpp"
#include "hdf5_io.hpp"
#include "misc.hpp"
#include "multicomp.hpp"

extern "C" void MeshRead(int *, int *, double *, char *);

template<typename T>
string ZeroPadNumber(T num){
    ostringstream ss;
    ss << setw( 5 ) << setfill( '0' ) << (int)num;
    return ss.str();
}

void scale_pos(Vec3d *pos, double R, int N){
  for(int i = 0; i<N; i++) pos[i] = pos[i]*R;
}

bool isPlaner(Vec3d *pos, int Np){
    for(int i=0;i<Np;i++){if(pos[i].z!=0){return false;}}
    return true;
}

pair<double, double> get_box_dim(Vec3d *Pos, MESH_p mesh){
    double xmin=1e7,xmax=-1e7,ymin=1e7,ymax=-1e7;
    for (int i = 0; i < mesh.N; ++i){
        if (Pos[i].x<xmin)  xmin=Pos[i].x;
        if (Pos[i].x>xmax)  xmax=Pos[i].x;
        if (Pos[i].y<ymin)  ymin=Pos[i].y;
        if (Pos[i].y>ymax)  ymax=Pos[i].y;
    }
    double xlen = xmax-xmin;
    double ylen = ymax-ymin;
    return {xlen, ylen};
}

double start_simulation(Vec3d *Pos, MESH_p &mesh, McP mcobj, STE &stretchobj, 
    string outfolder, double radius, int &residx){

    double Pole_zcoord;
    double ave_bond_len;
    int tdumpskip, titer;
    string resfile;

    hdf5_io_read_double( (double *)Pos,  outfolder+"/input.h5", "pos" );
    scale_pos(Pos, radius, mesh.N);
    hdf5_io_read_mesh((int *) mesh.numnbr, (int *) mesh.node_nbr_list, 
        outfolder+"/input.h5");

    if (isPlaner(Pos,mesh.N)){
        mesh.topology="flat";
        mesh.edge = get_nstart(mesh.N, 1);
        // To make sure that the edges do not coincide after pbc.
        mesh.boxlen=get_box_dim(Pos, mesh).first*(1+1/sqrt(mesh.N));
    }
    else{
        mesh.topology="sphere";
        mesh.edge = 0;
        mesh.boxlen=0;
        mesh.bdry_type=2;   // Sphere always has a pbc.
    }

    ave_bond_len = stretchobj.init_eval_lij_t0(Pos, mesh, mcobj.isfluid(),
        mesh.boxlen);
    if(!mcobj.isrestart()) { residx = 0; }
    else{
        fstream restartfile(outfolder+"/restartindex.txt", ios::in);
        restartfile >> titer >>  tdumpskip;
        restartfile.close();
        residx = titer;
        // if(area_para.do_area)init_area_t0(Pos,mesh,mbrane_para,area_para);
        // init_spcurv(spcurv_para, Pos, mbrane_para.N);
        // if(stick_para.do_stick)
            // identify_attractive_part(Pos, stick_para.is_attractive, stick_para.theta, mbrane_para.N);
        // max(&mesh.nPole,&Pole_zcoord,Pos,mbrane_para.N);
        // min(&mesh.sPole,&Pole_zcoord,Pos,mbrane_para.N);
        // identify_attractive_part(Pos, stick_para.is_attractive, stick_para.theta, mbrane_para.N);
        resfile=outfolder+"/snap_"+ZeroPadNumber(titer/tdumpskip)+".h5";
        hdf5_io_read_double( (double *)Pos,  resfile, "pos");
        hdf5_io_read_mesh((int *) mesh.numnbr,
                (int *) mesh.node_nbr_list, resfile);
    }
    // init_activity(act_para, mbrane_para.N);
    return ave_bond_len;
}
/*----------------------------------------------------------*/
void diag_wHeader(BE bendobj, STE steobj, std::fstream &fid ){
  std::string log_headers = "#iter acceptedmoves bend_e stretch_e ";
    // if(stick_para.do_stick){log_headers+="stick_e ";}
    // if(afm_para.do_afm){log_headers+="afm_e ";}
    // if (spring_para.do_spring){log_headers+="spring_e ";}
    if(steobj.dopressure()) {log_headers+=" Pressure_e ";}
    if(steobj.dovol()) {log_headers+=" Volume_e ";}
    // if(area_para.do_area) {log_headers+="area_e ";}
    // if (vol_para.is_pressurized){log_headers+="pressure_e ";}
    // if (afm_para.do_afm){log_headers+="afm_fx, afm_fy afm_fz ";}
    // if (spring_para.do_spring){log_headers+="spr_north.z spr_south.z ";}
    // /* log_headers+="volume nPole_z sPole_z hrms"; */
    log_headers+="total_e  volume  area";
    fid << log_headers << endl;
}
/*----------------------------------------------------------*/
int main(int argc, char *argv[]){
    int mpi_err, mpi_rank, residx;
    mpi_err = MPI_Init(0x0, 0x0);
    mpi_err =  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    pid_t pid = getpid();
    std::string fname;
    uint32_t seed_v;
    int iter, start, num_moves, num_bond_change, recaliter, num_exchange;
    double av_bond_len, Etot, radius;
    BE bendobj;
    STE stretchobj;
    MulCom lipidobj;
    McP mcobj(bendobj, stretchobj, lipidobj);
    Vec3d *Pos;
    MESH_p mesh;
    string outfolder, para_file, outfile, filename;
    char tmp_fname[128];

    start=0;
    outfolder = ZeroPadNumber(mpi_rank+start)+"/";
    fstream fileptr(outfolder+"/mc_log", ios::app);
    // Check if the file opened successfully
    //
    seed_v = (uint32_t) (mpi_rank + time(0));
    // seed_v=42
    RandomGenerator::init(seed_v);
    //
    para_file = outfolder+"/para_file.in";
    sprintf(tmp_fname, "%s", para_file.c_str());

    MeshRead(&mesh.bdry_type, &mesh.nghst, &radius, tmp_fname);
    mesh.N = (int)hdf5_io_get_Np(outfolder+"/input.h5", "pos")/3;
    Pos = (Vec3d *)calloc(mesh.N, sizeof(Vec3d));
    mesh.numnbr = (int *)calloc(mesh.N, sizeof(int));
    mesh.node_nbr_list = (int *)calloc(mesh.nghst*mesh.N, sizeof(int));

    mcobj.initMC(mesh.N, outfolder);
    bendobj.initBE(mesh.N, outfolder);
    stretchobj.initSTE(mesh.N, outfolder);
    mesh.av_bond_len = start_simulation(Pos, mesh, mcobj, stretchobj, outfolder, 
                radius, residx);
    lipidobj.initMulCom(mesh.N, mesh.av_bond_len, bendobj, outfolder);

    // How often do you want to compute the total energy?
    if (mcobj.isfluid()) recaliter=mcobj.fluidizeevery();
    else recaliter=10;
    
    clock_t timer;
    Etot = mcobj.evalEnergy(Pos, mesh);
    mcobj.write_energy(fileptr, iter);

    fstream outfile_terminal(outfolder+"/terminal.out", ios::app);
    outfile_terminal << "# The seed value is " << seed_v << endl;
    if(!mcobj.isrestart()) diag_wHeader(bendobj, stretchobj, fileptr);
    if(mcobj.isrestart()) fileptr << "# Restart index " << residx << endl;

    for(iter=residx; iter < mcobj.totaliter(); iter++){
        if(iter%mcobj.dumpskip() == 0){
            outfile=outfolder+"/snap_"+ZeroPadNumber(iter/mcobj.dumpskip())+".h5";
            hdf5_io_write((double*) Pos, 3*mesh.N, outfile, "pos");
            hdf5_io_write_mesh(mesh.numnbr, mesh.node_nbr_list,
                                mesh.N, mesh.nghst, outfile);
            hdf5_io_write(lipidobj.lipA.data(), mesh.N, outfile, "lip");
            fstream restartfile(outfolder+"/restartindex.txt", ios::out);
            restartfile << iter << " " << mcobj.dumpskip() << endl;
            restartfile.close();
        }
        num_moves = mcobj.monte_carlo_3d(Pos, mesh);
        num_exchange = mcobj.monte_carlo_lipid(Pos, mesh);
        if(!(iter % recaliter)){
            Etot = mcobj.evalEnergy(Pos, mesh);
            cout << "iter = " << iter << 
            "; Accepted Moves = " << (double)num_moves*100/mcobj.onemciter() 
            << " %;"
            "; Exchanged Moves = " << (double)num_exchange * 100 / mcobj.onemciter()
            << " %;"
            << " totalener = " << Etot << "; volume = " << mcobj.getvolume() << endl;
        }
        if (mcobj.isfluid() && !(iter % mcobj.fluidizeevery())) {
            num_bond_change = mcobj.monte_carlo_fluid(Pos, mesh);
            cout << "fluid stats " << num_bond_change << " bonds flipped" << endl;
        }
        mcobj.write_energy(fileptr, iter);
    }
    outfile_terminal << "Total time taken = " << (clock()-timer)/CLOCKS_PER_SEC << "s" << endl;
    outfile_terminal.close();
    fileptr.close();
    free(Pos);
    free(mesh.node_nbr_list);
    MPI_Barrier(MPI_COMM_WORLD);
    mpi_err = MPI_Finalize();
    return 0;
}