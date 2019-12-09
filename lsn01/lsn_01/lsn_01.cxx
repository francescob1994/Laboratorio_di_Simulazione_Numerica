/*#######################################
#                                       #
#  Leboratorio di Simulazione Numerica  #
#                                       #
#             Lezione 01                #
#                                       #
#          Francesco Bollati            #
#                                       #
########################################*/

#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include"random.h"

using namespace std;

void DataBlock(double[], double[], double[], int ,int);
double error(const double*, const double*, int );
void PrintData(char*, double[], double[], int, int ); 
double Chi2Test(double[], double[], int, int);

double ExpDist(double);
double CLDist(double);  
double DiceDist(double);
void CreateSample(Random, double (*)(double), char* );


              
int main (int argc, char *argv[]){


//************************************ ESERCIZIO 01.1 ************************************// 

   int M = 100000;	 //total number of throws
   int N = 100;          //number of blocks
   int L = int(M/N);     //number of throws in each block

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


// 1)

   //sample of r ~ U(0,1)

   double r[M];
   for(int i=0; i<M; i++) r[i] = rnd.Rannyu();   

   //Computation of blocks' means, progressive means and mean errors  
  
   double ave_prog[N], err_prog[N];
   DataBlock(r, ave_prog, err_prog, N, L);

   //print progressive means and errors on data1.dat

   char* file = "data1.dat";
   PrintData(file, ave_prog, err_prog, N, L);    


// 2)

   //Sample of sigma^2
   
   double sigma2[M];
   for(int i=0; i<M; i++)  sigma2[i] = pow( r[i] - 0.5 , 2 );

   //Computation of blocks' means, progressive means and mean errors  
  
   DataBlock(sigma2, ave_prog, err_prog, N, L);

   //print progressive means and errors on data2.dat

   file = "data2.dat";
   PrintData(file, ave_prog, err_prog, N, L);  


// 3)

   //partition interval [0,1]

   int Nbins = 100;
   double part[Nbins+1];            
   part[0] = 0;
   for(int i=1; i<Nbins+1; i++){
      part[i] = part[i-1] + 1./100;
   }

   //Chi2Test

   double Chi2;
   int n = 10000;       //dim sample
   double sample[n];
   file = "Chi2.dat";
   ofstream chi;
   chi.open(file);

   for(int i=0; i<100; i++){
          
         for(int j=0; j<n; j++) sample[j] = rnd.Rannyu();        
 
         Chi2 = Chi2Test(sample, part, Nbins, n);
         chi << i+1 << "   " << Chi2 << endl;
   }

   chi.close();         




//************************************ ESERCIZIO 01.2 ************************************//


   file = "Dice.dat";
   CreateSample( rnd, DiceDist, file );

   file = "Exp.dat";
   CreateSample( rnd, ExpDist, file );

   file = "CL.dat";
   CreateSample( rnd, CLDist, file );



//************************************ ESERCIZIO 01.3 ************************************//


   int l = 10;            //square's side length
   double ll = 0.5;       //needle length
   int d = 1;             //distance between horizontal lines
   int Mtot = 800000;     //number of throws
   double ya,yb,vx,vy;
   int Nh = 0;            //number of needles landing on a line
   double sin_theta;      //sin of the angle between the needle and the positive semiaxis of x
   double pi_sample[Mtot];
   int t;                 //first time a needle hit a line

   for(int i=0; i<Mtot; i++){

	//throw extremum A of the needle between line 1 and line l+1
	ya = rnd.Rannyu(1,l+1);

        //generate theta from uniform distribution in [0,2*pi] 
	vx = 1;
	vy = 1;
	while( vx*vx + vy*vy > 1){
		vx = rnd.Rannyu(-1,1);
		vy = rnd.Rannyu(-1,1);
	}
        sin_theta = vy/sqrt(vx*vx + vy*vy);
	yb = ya + ll*sin_theta;

	//check if needle cuts a line
	if( int(ya)!=int(yb) ) Nh++;       
        if(Nh>0){ 
           pi_sample[i] = (2*ll*(i+1))/(Nh*d);
           if(pi_sample[i-1]==0) t = i;
        }else{
           pi_sample[i] = 0;
        }  
   }

   Mtot = Mtot-t;
   L = int(Mtot/N);
   DataBlock(pi_sample, ave_prog, err_prog, N, L);

   file = "Pi.dat";
   PrintData(file, ave_prog, err_prog, N, L);  


   rnd.SaveSeed();
   return 0;
}


//FUNCTIONS

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

double Chi2Test(double sample[], double part[], int Nbins, int n){

   double chi2=0;
   double Ei = n/double(Nbins);
   for(int i=0; i<Nbins; i++){
      double Oi = 0;
      for(int j=0; j<n; j++){
         if(sample[j]>= part[i] && sample[j]<part[i+1]) Oi++;
      }
      chi2 += pow((Oi-Ei),2)/Ei;
   }
   return chi2;

}


double ExpDist(double u){  return -log(1.-u);  }
double CLDist(double u){  return tan( M_PI*(u-0.5) );  }
double DiceDist(double u){
   
    if( u<1./6 ){
       return 1;
    }else if( u>=1./6 && u<1./3 ){
       return 2;
    }else if( u>=1./3 && u<1./2 ){
       return 3;
    }else if( u>=1./2 && u<2./3 ){
       return 4;
    }else if( u>=2./3 && u<5./6 ){
       return 5;
    }else{
       return 6;
    }

}

void CreateSample(Random rnd, double (*dist)(double), char* file){

   ofstream write;
   write.open(file);

   int N[4] ={1, 2, 10, 100};
   int n = 10000;             //dim sample SN
   double sample[n][4];

   for(int i=0; i<4; i++){
      for(int j=0; j < n ; j++){
         sample[j][i] = 0;
	 for(int k=0; k<N[i]; k++){
            sample[j][i] += dist(rnd.Rannyu());
         }
         sample[j][i]/=N[i];
      }
   }
   for(int i=0; i<n; i++){
      write << sample[i][0] << "   " << sample[i][1] << "   " << sample[i][2] << "   " << sample[i][3] << endl;
   }

   write.close();

}
























