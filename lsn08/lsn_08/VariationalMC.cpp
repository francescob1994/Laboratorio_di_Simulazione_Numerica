/*#######################################
#                                       #
#  Leboratorio di Simulazione Numerica  #
#                                       #
#             Lezione 08                #
#                                       #
#          Francesco Bollati            #
#                                       #
########################################*/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "VariationalMC.h"

using namespace std;

double epsilon = 0.13;

int main()
{ 

  Input(); //Inizialization

  //FindEquilibrationTime();
  //Correlations();

  Equilibration();

  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
      Reset(iblk);   //Reset block averages
      for(int istep=1; istep <= nstep; ++istep)
      {
         Move();
         Measure();
         if( istep % 10 == 0 ) Print_x();         
         Accumulate(); //Update block averages
       }
       Averages(iblk);   //Print results for current block
  }

  PrintFinalH();

  return 0;
}




void Input(void)
{
  ifstream ReadInput;

  cout<<" Quantum particle 1D                   " <<endl<<endl;
  cout<<" Variational MonteCarlo Simulation     " <<endl<<endl;
  cout<<" V(x) = x^2 + 5*x^4/2                  " <<endl<<endl;
  cout<<" Psi_Trial = Gaussian( -mu, sigma ) + Gaussian( +mu, sigma )              " <<endl<<endl;
  cout<<" The program uses a0=1, hbar=1 and m=1" <<endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");
  
  ReadInput >> mu;
  ReadInput >> sigma;
  cout << " Parameters of the trial wave function:    Mu = " << mu << "    Sigma = " << sigma << endl << endl;

  norm = 2.*sigma*sqrt(pi)*( 1 + exp(-(mu*mu)/(sigma*sigma)) ); 
    
  ReadInput >> delta;
  //delta = sigma + 0.5;

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> teq; 

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Equilibration time = "<< teq << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;

  ReadInput >> x0; 
  x = x0;

  cout << "Starting point of the system x0 = " << x0 << endl << endl;


  ReadInput.close();
}


double e1(double x){ return exp(-(x+mu)*(x+mu)/(2*sigma*sigma)); }
double e2(double x){ return exp(-(x-mu)*(x-mu)/(2*sigma*sigma)); }
double psi(double x){ return e1(x) + e2(x); }
double V(double x){ return pow(x,4)-5.*x*x/2.; }
double D2psi(double x){ return -psi(x)/(sigma*sigma) + ( e1(x) * pow((x+mu),2) + e2(x) * pow((x-mu),2) )/pow(sigma,4); }
double Hpsi(double x){ return -0.5*D2psi(x) + V(x)*psi(x); }
double integrand(double x){ return Hpsi(x)/psi(x); }
double rho(double x){ return psi(x)*psi(x)/norm; }

void Move(void)
{ 
   double p, xold, xnew;
   double dx = delta*(2*rnd.Rannyu() -1); //random number in (-delta, delta)
   xold = x;
   xnew = x + dx;

   p = rho(xnew)/rho(xold);

   if(p >= rnd.Rannyu() ){
        x = xnew;
        accepted++;
   }
   attempted++;            
};


void Measure(void)
{
  oss = integrand(x);
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       glob_av = 0;
       glob_av2 = 0;
   }

   blk_av = 0;
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}



void Accumulate(void) //Update block averages
{

   blk_av = blk_av + oss;
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream H;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    H.open("output.H.0",ios::app);

    stima_H = blk_av/blk_norm;
    glob_av += stima_H;
    glob_av2 += stima_H*stima_H;
    err_H = Error(glob_av, glob_av2, iblk);

    if(iblk%100 == 0){
       H << setw(wd) << iblk <<  setw(wd) << stima_H << setw(wd) << glob_av/(double)iblk << setw(wd) << err_H << endl;
    }
    cout << "----------------------------" << endl << endl;

    H.close();
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}


void FindEquilibrationTime(void){
 
    ofstream eq;
    eq.open("equilibration.dat");
    double o;

    for(int t=0; t<600; t++){

        o = 0;
	Measure();
        o = oss;
        eq<<o<<"   "<<x<<endl; 

        Move();
    }
    eq.close();

};


void Equilibration(void){

    for(int t=0; t<teq; t++) Move();
   
    //cout<<"System Equilibrated"<<endl<<endl;

}


void Correlations(void)  //per stabilire il valore di nsteps
{

   ofstream cor;
   cor.open("correlation.dat");
   
   x = x0;
   int Omega = 1000;
   double X; 

   int tauf = 200;
   double cov[tauf];
   for(int t=0; t<tauf; t++) cov[t]=0;

   for(int w = 1; w <= Omega; w++){ //ciclo sulle traiettorie

        Equilibration();
        X = x;

   	for(int t = 1; t <= tauf; t++){

                Move();
                cov[t-1] += X*x/Omega;
                
                if(w == Omega) cor<<t<<"   "<<cov[t-1]<<endl; 
        }        
    }
    cor.close();

}    


void PrintFinalH(void)
{

   ofstream Hfin;
   Hfin.open("finalH.dat", ios::app);
   //const int wd = 12;
  
   Hfin << mu << "   " << sigma << "   " << glob_av/double(nblk) << "   " << err_H << endl;

}

void Print_x(void)
{

   ofstream xx;
   xx.open("x.dat", ios::app);
   xx << x << endl;
   xx.close();

}
