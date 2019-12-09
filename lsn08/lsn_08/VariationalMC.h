#ifndef __variationalMC_
#define __variationalMC_

#include<cmath>

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
double oss;     //value of integrand = Hpsi/psi in x

// averages
double blk_av, blk_norm, accepted, attempted;
double glob_av, glob_av2;
double stima_H, err_H;

//pigreco
const double pi=3.1415927;

//psi_trial paramenters
double mu, sigma;
double norm;     //normalization probability density

//state of the system
double x0, x;

//simulation
int nstep, nblk;
double delta;
int teq; 

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void Measure(void);
double Error(double,double,int);

void FindEquilibrationTime(void);
void Equilibration(void);
void Correlations(void);
void PrintFinalH(void);
void Print_x(void);

double psi(double);
double V(double);
double D2psi(double);
double Hpsi(double);
double integrand(double);
double rho(double);

#endif
