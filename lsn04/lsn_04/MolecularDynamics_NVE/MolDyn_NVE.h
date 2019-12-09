/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=4;
int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres; //aggiunto pressione

//g(r)
int nbins;
double walker[100], blk_gofr[100], glob_gofr[100], glob_gofr2[100];
double bin_size, stima_gofr, err_gofr;

// averages
double acc,att;

//configuration
const int m_part=108;  
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;
bool old;

//blocchi //
int nblocks, StepsPerBlock, n_mis, iblk;
double v_blk, k_blk, t_blk, e_blk, p_blk;
double v_prog, k_prog, t_prog, e_prog, p_prog;
double v_prog2, k_prog2, t_prog2, e_prog2, p_prog2;
double v_err, k_err, t_err, e_err, p_err; 

//pigreco
const double pi=3.1415927;
  
//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfOld(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);

double Error(double, double, int);
void Reset(int);
void Accumulate(void);
void MeasureAve(int);

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
