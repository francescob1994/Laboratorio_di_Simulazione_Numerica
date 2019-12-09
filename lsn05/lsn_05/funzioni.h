#ifndef __funzioni_h__
#define __funzioni_h__

#include<iostream>
#include<cmath> 
#include"random.h"

using namespace std;


//STRUTTURA DATI

struct punto {

	double x;
        double y;
        double z;

        //spherical coordinates
        double r(void){ return sqrt(x*x + y*y + z*z); }
        double theta(void){ 
        	double rho = sqrt(x*x + y*y);
                if(rho > 0){
                	if(z > 0) return atan(rho/z); 
                	if(z == 0) return M_PI/2.;
                	if(z < 0) return M_PI + atan(rho/z);
                }
                if(rho == 0){
                	if(z > 0) return 0;
                        if(z < 0) return M_PI;
               }

        }
        double phi(void){
                double rho = sqrt(x*x + y*y);
                if(y > 0) return acos(x/rho);
                if(y < 0) return 2*M_PI - acos(x/rho);

        }
        
        //scalar product
        double dot( const punto & P ){ return x*P.x + y*P.y + z*P.z; }
       
        //overloading logic operator ==
        bool operator==(const punto & P){
        	if( x == P.x && y == P.y && z == P.z ) return true;
                else return false;
        } 
};


//PROBABILITY DENSITIES

//(psi*)=psi*a0^(3/2), (r*)=r/a0.
double psi2_100(punto);
double psi2_210(punto);

//METROPOLIS ALGORITHMS

void Metropolis1(punto&, double, Random&, double (*rho)(punto)); 
void Metropolis2(punto&, double, Random&, double (*rho)(punto));

//DATABLOCK

void DataBlock(double[], double[], double[], int ,int);
double error(const double*, const double*, int );
void PrintData(char*, double[], double[], int, int ); 

//MONTECARLO SIMULATION

void MonteCarlo(int, int, void (*Metropolis)(punto&, double, Random&, double (*rho)(punto)), double (*rho)(punto), int, punto, Random&, char*, double, double[], double[], double[]);


#endif




