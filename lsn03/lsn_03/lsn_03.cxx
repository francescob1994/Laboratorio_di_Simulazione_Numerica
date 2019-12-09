/*#######################################
#                                       #
#  Laboratorio di Simulazione Numerica  #
#                                       #
#             Lezione 03                #
#                                       #
#          Francesco Bollati            #
#                                       #
########################################*/

#include<iostream>
#include<string>
#include"random.h"
#include<cmath>
#include<fstream>

using namespace std;


void DataBlock(double[], double[], double[], int ,int);
double error(const double*, const double*, int );
void PrintData(char*, double[], double[], int, int ); 

//functions for Black-Sholes analytic solution 
double N(double);
void black_scholes(double, double, double, double, double, double&, double&);


int main( int argc, char *argv[] ){

  
   //Random object
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;



  
   int N = 100;              //number blocks
   int Nstep = 101;          //number time steps
   int M = 10000;            //dim dataset for each t
   int L = int(M/N);
   double S0 = 100;          //asset price(t=0)                          
   double T = 1;             //delivery time
   double dt = T/(Nstep -1); //time interval between consecutive time steps
   double K = 100;           //strike price
   double r = 0.1;           //risk-free interest rate
   double sigma = 0.25;      //volatility


   // DIRECT SAMPLING OF S(T)
   
   double ST = 0;

   double C[M];              //Call sample
   double errC[N];           //Call block errors
   double meanC[N];          //Call progressive means

   double P[M];              //Put sample
   double errP[N];           //Put block errors
   double meanP[N];          //Put progressive means

   for(int w = 0; w<M; w++){

   	ST = S0*exp( (r-0.5*sigma*sigma)*T + sigma*rnd.Gauss(0,T) );

        C[w] = exp(-r*T)*max(0., ST - K );
        P[w] = exp(-r*T)*max(0., K - ST );
   } 
   
   DataBlock(C, meanC, errC, N, L);
   DataBlock(P, meanP, errP, N, L);

   //print data

   char* fileC = "Call_direct.dat";
   char* fileP = "Put_direct.dat";
   PrintData(fileC, meanC, errC, N, L);
   PrintData(fileP, meanP, errP, N, L);

 

   // SAMPLING PATH OF S(t)

   //allocate memory for sample of trajectories
   double **S = new double *[M];                           
   for(int w=0; w<M; w++){
   	S[w] = new double[Nstep];	
   }

   //each trajectory begins with S0
   for(int w=0; w<M; w++){  
   	S[w][0] = S0;       
   }

   //generate trajectories

   for(int w=0; w<M; w++){

   	for(int i=0; i<Nstep-1; i++){
           S[w][i+1] = S[w][i]*exp( (r-0.5*sigma*sigma)*dt + sigma*rnd.Gauss(0,sqrt(dt)) );                     
        }

        C[w] = exp(-r*T)*max(0., S[w][100] - K );
        P[w] = exp(-r*T)*max(0., K - S[w][100] );
   }  


   DataBlock(C, meanC, errC, N, L);
   DataBlock(P, meanP, errP, N, L);

   //print some trajectories

   ofstream traj;
   traj.open("paths.dat");
   
   for(int t=0; t<Nstep; t++){
         traj << t*dt << "   ";
      for(int w=0; w<M; w+=100){
         traj << S[w][t] << "   "; 
      }
      traj<<endl;
   }

   traj.close();

   //print data

   fileC = "Call_path.dat";
   fileP = "Put_path.dat";
   PrintData(fileC, meanC, errC, N, L);
   PrintData(fileP, meanP, errP, N, L);
   

    
   //ANALYTICAL SOLUTION
   double Call, Put;
   black_scholes(S0, K, T, r, sigma, Call, Put);
   cout<<"BLACK-SCHOLES:   Call = "<<Call<<",   Put = "<<Put<<endl;

 
   delete[] S;
   rnd.SaveSeed();
   return 0;
}


//FUNCTIONS



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


double error(const double* AV, const double* AV2, int n){

   if(n==0){ return 0;
   }else{
	return sqrt( (AV2[n] - AV[n]*AV[n])/n );
   }
};




double N( double x){
    return 0.5 * (1. + erf(x / sqrt(2.)));
};

void black_scholes(double S0, double K, double T, double r, double sigma, double & C, double & P){
    double d1 = 1./(sigma * sqrt(T)) * (log(S0 / K) + (r + (sigma*sigma) / 2.) * T);
    double d2 = d1 - sigma * sqrt(T);
    C = S0 * N(d1) - K * exp(-r * T) * N(d2);
    P = S0 *(N(d1) - 1.) - K * exp(-r * T) * (N(d2)-1.);
};


