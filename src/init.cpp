#include "global.h"
#include "subroutine.h"
/* std::mt19937 rng2; */

extern "C"  void Membrane_listread(int *, double *,  double *, 
            double *, int *, char *);

extern "C"  void Area_listread(bool *, double *, double *, char *);

extern "C"  void Stick_listread(bool *, double *, double *, double *, 
            double *,  char *);

extern "C"  void  MC_listread(char *, double *, double *, bool *,
            int *, int *, char *);

extern "C"  void  Activity_listread(char *, double *, double *, char *);

extern "C"  void  Afm_listread(bool *, double *, double *, double *, double *,
             char *);

extern "C"  void  Spring_listread(bool *, double *, double *, double *, char *);

extern "C"  void  Fluid_listread(bool *, int * , int *, double *, char *);
extern "C"  void   Volume_listread(bool *, bool *, double *, double*, char *); 

/* void init_rng2(uint32_t seed_val) { */

  ///  @brief Generates random number
  /* rng2.seed(seed_val); */
/* } */



void init_system_random_pos(Vec2d *Pos,  double len, 
        int N, char *metric, int bdry_condt ){

    /// @brief Initializes the points on surface of sphere or flat plane
    ///  @param Pos array containing co-ordinates of all the particles
    ///  @param metric Topology of the surface "cart" for flat plane "sph" for
    /// sphere
    ///  @param len length of the domain;
    ///  @param N number of points;

    bool is_sph, is_cart;
    double dl;

    // this should be calculated or passed parameters
    int n_ghost;
    // remove it once debugged 

    n_ghost = (int) sqrt(N);
    dl = (len/n_ghost);

    /* std::uniform_real_distribution<> rand_real(dl, len-dl); */

    
    is_sph = false;
    is_cart = false;
    if(strcmp(metric, "sph") == 0){
        is_sph = true;
    }
    if(strcmp(metric, "cart") == 0){
        is_cart = true;
    }

    if(!is_cart && !is_sph){
        fprintf(stderr, "Unknown metric, input should be: \n");
        fprintf(stderr, "a) cart for a 2d cartesian monte carlo \n");
        fprintf(stderr, "b) sph  for a monte carlo on the surface of a sphere\n");
        exit(0);
    }

     //::generateUniform  rngu; //(0.0, 1.0)
     
    if(is_cart){
        switch (bdry_condt) {
            case 0:
            // This is a channel ; read bdry_condt as 0 
            // n_ghost has to be even for logic to work;
            n_ghost = 2*n_ghost;
            for(int i=0; i<n_ghost/2; i++){
                    Pos[i].x = i*dl;
                    Pos[i].y = 0.0; 
            }
            for(int i=n_ghost/2; i<n_ghost; i++){
                Pos[i].x = (i - n_ghost/2 + 0.5)*dl;
                /* Pos[i].x = (i - n_ghost/2)*(2*(len + 0.5)/n_ghost);  bp*/
                Pos[i].y = len; 
            }
            for(int i=n_ghost; i<N; i++){
                Pos[i].x = RandomGenerator::generateUniform(dl, len-dl);//rand_real(rng2);
                Pos[i].y = RandomGenerator::generateUniform(dl, len-dl); //rand_real(rng2);
            }
            break;
      case 1:
            // This is a frame ;
            // n_ghost has to be even for logic to work;
            n_ghost = 4*n_ghost;
            for(int i=0; i<n_ghost/4; i++){
                Pos[i].x = (i+0.5)*dl;
                Pos[i].y = 0.0; 
            }
            for(int i=n_ghost/4; i<n_ghost/2; i++){
                Pos[i].x = (i - n_ghost/4 + 0.5)*dl;
                Pos[i].y = len; 
            }
            for(int i=n_ghost/2; i<3*n_ghost/4; i++){
                Pos[i].x = 0.0;
                Pos[i].y = (i - n_ghost/2 + 0.5)*dl; 
            }

            for(int i=3*n_ghost/4; i<n_ghost; i++){
                Pos[i].x = len;
                Pos[i].y = (i +0.5 - 3*n_ghost/4)*dl; 
            }
            for(int i=n_ghost; i<N; i++){
                Pos[i].x = RandomGenerator::generateUniform(dl, len-dl);//rand_real(rng2);
                Pos[i].y = RandomGenerator::generateUniform(dl, len-dl); //rand_real(rng2);
            }

            break;

      default:
            for(int i=0; i<N; i++){
                Pos[i].x = RandomGenerator::generateUniform(dl, len-dl);//rand_real(rng2);
                Pos[i].y = RandomGenerator::generateUniform(dl, len-dl); //rand_real(rng2);
            }
        }
    }
    if(is_sph){
        Pos[0].x = 0;
        Pos[0].y = 0;
        Pos[1].x = pi;
        Pos[1].y = 0;
        for(int i=2; i<N; i++){
            Pos[i].x = acos(2*drand48() - 1); 
            Pos[i].y = 2*pi*drand48();
        }
        Pos[2].x = acos(2*drand48() - 1); 
        Pos[2].y = 2*pi*drand48();

    }
    /* for(int i=0; i<N; i++){ */
    /*     printf("%lf %lf \n", Pos[i].x); */
    /*     Pos[i].x = drand48()*len; */
    /*     Pos[i].y = drand48()*len; */
    /* } */

}

void init_eval_lij_t0(Vec3d *Pos, MESH_p mesh, double *lij_t0,
         MBRANE_p *para, SPRING_p *spring, bool is_fluid){
    /// @brief evaluates distance between neighbouring points and stores in lij_t0
    ///  @param Pos array containing co-ordinates of all the particles
   /// @param lij_t0 initial distance between points of membrane
    ///  @param mesh mesh related parameters -- connections and neighbours
    /// information; 
    ///  @param para membrane related parameters 
    Vec3d dr;
    int i,j,k;
    int num_nbr, cm_idx, npairs;
    double sum_lij=0;
    double r0;
    npairs = 0;
    for(i = 0; i < para->N; i++){
        num_nbr = mesh.numnbr[i];
        cm_idx = mesh.nghst * i;
        for(k = cm_idx; k < cm_idx + num_nbr; k++) {
            j = mesh.node_nbr_list[k];
            /* dr = diff_pbc(Pos[i] , Pos[j], para->len); */
            dr = Pos[j] - Pos[i];
            lij_t0[k] = sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
            sum_lij += sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
            npairs++;
            /* printf("%g %g %g %g %g \n", Pos[i].x, Pos[j].x, Pos[i].y, Pos[j].y, lij_t0[k]); */
        }
    }
    para->av_bond_len = sum_lij/npairs;
    r0=para->av_bond_len;
    spring->constant=para->coef_bend/(r0*r0);
    if(is_fluid){
        for(i = 0; i < mesh.nghst*para->N; i++){
            lij_t0[i] = para->av_bond_len;
        }
    }
}

void init_read_parameters(MBRANE_p *mbrane_para, MC_p *mc_para, AREA_p *area_para, 
        FLUID_p *fld_para, VOL_p *vol_para, STICK_p *stick_para, AFM_p *afm_para, 
        ACTIVE_p *act_para, SPRING_p *spring_para, LIPID_p *lipid_para,
        string para_file){
   /// @brief read parameters from para_file 
    ///  @param mesh mesh related parameters -- connections and neighbours
    /// information; 
    ///  @param mbrane membrane related parameters
    ///  @param afm AFM related parameters
    ///  @param mcpara monte-carlo related parameters
    //
    //
    char temp_algo[char_len];
    char which_act[char_len], tmp_fname[char_len];

    sprintf(tmp_fname, "%s", para_file.c_str() );
    Membrane_listread(&mbrane_para->N, &mbrane_para->coef_bend,
             &mbrane_para->sp_curv,  &mbrane_para->radius, 
            &mbrane_para->bdry_type, tmp_fname);

    bool do_stick;
    double pos_bot_wall;  // position of the bottom attractive wall
    double sigma, epsilon, theta; // sigma and epsilon for the bottom attractive wall
 
    Stick_listread(&stick_para->do_stick, &stick_para->pos_bot_wall, 
            &stick_para->sigma, &stick_para->epsilon, &stick_para->theta,
            tmp_fname);

    /* sprintf(tmp_fname, "%s", para_file.c_str() ); */
    MC_listread(temp_algo, &mc_para->dfac, &mc_para->kBT, &mc_para->is_restart,
            &mc_para->tot_mc_iter, &mc_para->dump_skip, tmp_fname);
    mc_para->algo=temp_algo;

    /* sprintf(tmp_fname, "%s", para_file.c_str() ); */
    Spring_listread(&spring_para->do_spring, &spring_para->constant,
            &spring_para->nPole_eq_z, &spring_para->sPole_eq_z, tmp_fname);

    /* sprintf(tmp_fname, "%s", para_file.c_str() ); */
    Afm_listread(&afm_para->do_afm, &afm_para->tip_rad, 
            &afm_para->tip_pos_z, &afm_para->sigma, &afm_para->epsilon,
             tmp_fname);

    /* sprintf(tmp_fname, "%s", para_file.c_str() ); */
    Activity_listread(which_act, &act_para->minA, &act_para->maxA, tmp_fname);
    act_para->act = which_act;

    sprintf(tmp_fname, "%s", para_file.c_str() );
    Fluid_listread(&fld_para->is_fluid, &fld_para->min_allowed_nbr,
            &fld_para->fluidize_every, &fld_para->fac_len_vertices, tmp_fname);

    /* sprintf(tmp_fname, "%s", para_file.c_str() ); */
    Volume_listread(&vol_para->do_volume, &vol_para->is_pressurized,
           &vol_para->coef_vol_expansion, &vol_para->pressure, tmp_fname);

    Area_listread(&area_para->is_tether, &area_para->YY,
            &area_para->Ka, tmp_fname);

    Lipid_listread(&lipid_para->islipid);
    mc_para->one_mc_iter = 2*mbrane_para->N;
    mc_para->delta = sqrt(8*pi/(2*mbrane_para->N-4));
}
//
//
void write_parameters(MBRANE_p mbrane, MC_p mc_para, AREA_p area_para, FLUID_p fld_para, 
        VOL_p vol_p, STICK_p stick_para, AFM_p afm_para,  ACTIVE_p act_para, 
        SPRING_p spring_para, string out_file){
 
    double FvK = area_para.YY*mbrane.radius*mbrane.radius/mbrane.coef_bend;
    ofstream out_;
    out_.open( out_file );
    out_<< "# =========== Model Parameters ==========" << endl
            << " Foppl von Karman number: FvK = " << FvK << endl
            << " Elasto-thermal number: ET = " << mc_para.kBT/mbrane.coef_bend*sqrt(FvK) << endl
            << " average bond length: r0 = " << mbrane.av_bond_len << endl;

    out_<< "# =========== Membrane Parameters ==========" << endl
            << " coef bend = " << mbrane.coef_bend << endl
            << " sp_curv " << mbrane.sp_curv << endl
            << " N " << mbrane.N << endl
            << " av_bond_len " << mbrane.av_bond_len << endl
            << " bdry_type " << mbrane.bdry_type << endl
            << " radius " << mbrane.radius << endl;

    out_<< "# =========== Monte Carlo Parameters ==========" << endl
            << " algo = " << mc_para.algo << endl
            << " dfac " << mc_para.dfac << endl
            << " dump_skip " << mc_para.dump_skip << endl
            << " kbt " << mc_para.kBT << endl
            << " delta " << mc_para.delta << endl
            << " is_restart " << mc_para.is_restart << endl
            << " tot_mc_iter " << mc_para.tot_mc_iter << endl
            << " one mc iter " << mc_para.one_mc_iter << endl;

    out_<< "# =========== Area Parameters ==========" << endl
            << " is tether = " << area_para.is_tether << endl
            << " YY " << area_para.YY << endl
            << " sigma " << area_para.Ka << endl;

    out_<< "# =========== Activity Parameters ==========" << endl
            << " which activity = " << act_para.act << endl
            << " minA = " << act_para.minA << endl
            << " maxA = " << act_para.maxA << endl;

    out_<< "# =========== Fluid Parameters ==========" << endl
            << " is fluid= " << fld_para.is_fluid << endl
            << " min_allowed_nbr = " << fld_para.min_allowed_nbr << endl
            << " fluid iter every " << fld_para.fluidize_every << endl
            << " factor_len_vertices = " << fld_para.fac_len_vertices << endl;

    out_<< "# =========== Volume Parameters ==========" << endl
            << " do volume= " << vol_p.do_volume << endl
            << " is pressurized = " << vol_p.is_pressurized << endl
            << " coef_vol_expansion " << vol_p.coef_vol_expansion << endl
            << " pressure  " << vol_p.pressure << endl;

    out_<< "# =========== Sticking Parameters ==========" << endl
            << " do stick " << stick_para.do_stick << endl
            << " pos_bot_wall  " << stick_para.pos_bot_wall << endl
            << " sigma " << stick_para.sigma << endl
            << " epsilon " << stick_para.epsilon << endl
            << " theta " << stick_para.theta << endl;

    out_<< "# =========== Spring Parameters ==========" << endl
            << " do spring " << spring_para.do_spring << endl
            << " constant " << spring_para.constant << endl
            << " nPole_eq_z " << spring_para.nPole_eq_z << endl
            << " sPole_eq_z " << spring_para.sPole_eq_z << endl;
    out_.close();
}
//
void init_KK_0(double *KK, AREA_p area_p, MESH_p mesh, int N){
    int i, k;
    int num_nbr, cm_idx;
    double kk_t = area_p.YY*sqrt(3)/2;
    int mpi_err, mpi_rank;
    // double rndg;
    /* std::cout << "Gaussian Random Number: " << gaussianNumber << std::endl; */
    mpi_err =  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    /* std::uniform_real_distribution<> rand_real(-kk_t/4.0, kk_t/4); */
    for(i = 0; i < N; i++){
        num_nbr = mesh.numnbr[i];
        cm_idx = mesh.nghst * i;
        for(k = cm_idx; k < cm_idx + num_nbr; k++) {
            // rndg = RandomGenerator::generateGaussian(0.0,1.0);
            KK[k] = kk_t;
            // + sqrt(kk_t)*rndg*0.5; //rand_real(rng2);
            /* printf("%d  %lf \n", mpi_rank, KK[k]); */
        }
    }
    /* exit(0); */
}
//
void init_activity(ACTIVE_p activity, int N){
    int i;
    if(activity.act == "random"){
        for(i=0;i<N;i++) activity.activity[i] = RandomGenerator::generateUniform(activity.minA,activity.maxA); 
    }
    if(activity.act == "constant"){
        for(i=0;i<N;i++) activity.activity[i] = activity.maxA;
    }
}
//
void init_lipidcomp(bool *lipA, int N){
    int i;
    for(i=0;i<N;i++) lipA[i] = RandomGenerator::intUniform(0,1);
}