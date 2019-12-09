#include<iostream>
#include<cmath>
#include<fstream>

#include"random.h"
#include"funzioni.h"


using namespace std;


//PROBABILITY DENSITIES

double psi2_100(punto R){ return exp(-2*R.r())/M_PI; };
double psi2_210(punto R){ return 1./(32.*M_PI)*pow(R.r(),2)*exp(-R.r())*pow(cos(R.theta()),2); };

//METROPOLIS ALGORITHMS

void Metropolis1( punto & R, double r, Random & rnd, double (*rho)(punto p) ){ 

	double theta = M_PI*rnd.Rannyu();
        double phi = 2*M_PI*rnd.Rannyu();
        punto P = R; 
        P.x += r*sin(theta)*cos(phi);
        P.y += r*sin(theta)*sin(phi);
        P.z += r*cos(theta);

        //acceptance
        if( (*rho)(P) > (*rho)(R) ){
        	R = P;
        }else{
        	if( rnd.Rannyu() < (*rho)(P) / (*rho)(R) ) R = P; 
        }
      
};

void Metropolis2(punto & R, double r, Random & rnd, double (*rho)(punto p)){ 


        punto P = R;
        P.x += rnd.Gauss(0,r);
        P.y += rnd.Gauss(0,r);
        P.z += rnd.Gauss(0,r);

        //acceptance
        if( (*rho)(P) > (*rho)(R) ){
        	R = P;
        }else{
        	if( rnd.Rannyu() < (*rho)(P) / (*rho)(R) ) R = P; 
        }

};

//DATABLOCK

double error(const double* AV, const double* AV2, int n){

   if(n==0){ return 0;
   }else{
	return sqrt( (AV2[n] - AV[n]*AV[n])/n );
   }
};


void DataBlock(double sample[], double ave_prog[], double err_prog[], int N ,int L){
   
   double ave[N], ave2[N], sum2_prog[N];

   for(int i=0; i<N; i++){           //cycle on blocks

   	double sum = 0;
	for(int j=0; j<L; j++){      //cycle on throws within a block	
		sum += sample[j+i*L];
	}
        ave[i] = sum/L;              //i-th block mean
        ave2[i] = ave[i]*ave[i];     //i-th block square mean
        
   }

   for(int i=0; i<N; i++){

        ave_prog[i] = 0;
        sum2_prog[i] = 0;

	for(int j=0; j<i+1; j++){          
		
        	ave_prog[i] += ave[j];     //sum of blocks' means up to i-th block
		sum2_prog[i] += ave2[j];  
        }
	ave_prog[i] /= (i+1);
	sum2_prog[i] /= (i+1);
	err_prog[i] = error(ave_prog,sum2_prog,i); 
   }
};

void PrintData(char* file, double av[], double err[], int N, int L){

   ofstream w;
   w.open(file);
   for(int i=0; i<N; i++){
      w << (i+1)*L << "   " << av[i] << "   " << err[i] <<endl;
   }
   w.close();

};

//MONTECARLO SIMULATION

void MonteCarlo(int Nsteps, int teq, void (*Metropolis)(punto & R, double r, Random & rnd, double (*rho)(punto p)), double (*rho)(punto p), int N_blocks, punto R0, Random & rnd, char* file, double r, double sample[], double mean[], double MeanError[]){

   ofstream points;
   points.open("points_210.dat");
   punto R = R0;
   for(int t=0; t<(Nsteps + teq + 1); t++){

      (*Metropolis)(R,r, rnd, (*rho)); 

      if( t > teq ){

          sample[t-teq -1] = R.r();  

          //print Markov chain's points 
          if( (t-teq-1)%200 == 0 ) points << R.x << "   " << R.y << "   " << R.z << endl;          
          

      }
   }
   points.close();
   int L = int(Nsteps/N_blocks);
   DataBlock( sample, mean, MeanError, N_blocks, L);
   PrintData(file, mean, MeanError, N_blocks, L);


};




