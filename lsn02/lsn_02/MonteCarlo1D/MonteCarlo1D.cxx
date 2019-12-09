/*#######################################
#                                       #
#  Leboratorio di Simulazione Numerica  #
#                                       #
#             Lezione 02                #
#                                       #
#          Francesco Bollati            #
#                                       #
########################################*/

#include<iostream>
#include<string>
#include<fstream>
#include<cmath>
#include"random.h"

using namespace std;


double g(double);                         //integrand function
double gRpositive(double);                //function equal to g in [0,1] and equal to zero in (1,+inf)

double d(double);                         //probability distribution, importance sampling 1
double f(double);                         //integrand function, importance sampling 1

double gauss(double, double, double);     //probability distribution, importance sampling 2
double f2(double);                        //integrand function, importance sampling 2

void DataBlock(double[], double[], double[], int ,int);
double error(const double*, const double*, int );
void PrintData(char*, double[], double[], int, int );


int main(int argc, char* argv[]){

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


   int M = 100000;        //total number of throws
   int N = 100;          //number of blocks
   int L = int(M/N);     //number of throws in each block

   double mean[N];
   double MeanError[N];   
   double integranda[M];

//************************************ CASO 1 ************************************//  

   //integrand sample

   for(int i=0; i<M; i++){
   	integranda[i] = g(rnd.Rannyu());
   }
   
   //Computation of blocks' means, progressive means and mean errors
 
   DataBlock(integranda, mean, MeanError, N, L); 

   //print progressive means and errors on data1.dat

   char* file = "data1.dat";
   PrintData(file, mean, MeanError, N, L); 
   
   cout<<"(NO I.S.) I = "<<mean[N-1]<<" +/- "<<MeanError[N-1]<<endl;


//************************************ CASO 2 ************************************//

// IMPORTANCE SAMPLING 1

   //integrand sample 

   int count = 0; 
   double u = 0;
   double r = 0;
   while( count!=M ){
   	u = rnd.Rannyu();
   	r = rnd.Rannyu();
        if(r < 1 - M_PI*M_PI*u*u/8. ){                      
		integranda[count] = f(u);
		count++;
	}
   }

   //Computation of blocks' means, progressive means and mean errors
 
   DataBlock(integranda, mean, MeanError, N, L); 

   //print progressive means and errors on data1.dat

   file = "data2.dat";
   PrintData(file, mean, MeanError, N, L); 

   cout<<"(I.S. 1) I = "<<mean[N-1]<<" +/- "<<MeanError[N-1]<<endl;



// IMPORTANCE SAMPLING 2

   //integrand sample 

   for(int i=0; i<M; i++){
      r = rnd.Gauss(0,0.35);
      if(r<0){ r = -r; }
      integranda[i] = f2(r);
   }

   //Computation of blocks' means, progressive means and mean errors
 
   DataBlock(integranda, mean, MeanError, N, L); 

   //print progressive means and errors on data1.dat

   file = "data3.dat";
   PrintData(file, mean, MeanError, N, L); 

   cout<<"(I.S. 2) I = "<<mean[N-1]<<" +/- "<<MeanError[N-1]<<endl;

 
   rnd.SaveSeed();
   return 0;
}


//FUNCTIONS



double g( double x ){ return M_PI*0.5*cos(M_PI*x/2); } 

double gRpositive( double x ){                         
   if(x<=1){ 
      return g(x);
   }else{
      return 0; 
   }
} 

double d( double x ){ return ( 24. - 3*M_PI*M_PI*x*x )/( 24.-M_PI*M_PI ); } 
double f( double x ){ return g(x)/d(x); }                                 

double gauss(double x, double mu, double sigma){ return 1./(sigma*sqrt(2*M_PI))*exp( -(x-mu)*(x-mu)/(2.*sigma*sigma) ); }
double f2( double x ){ return gRpositive(x)/(2*gauss(x,0,0.35)); }  

		      
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


double error(const double* AV, const double* AV2, int n){

   if(n==0){ return 0;
   }else{
	return sqrt( (AV2[n] - AV[n]*AV[n])/n );
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

