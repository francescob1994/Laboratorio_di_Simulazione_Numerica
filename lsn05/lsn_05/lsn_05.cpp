/*#######################################
#                                       #
#  Leboratorio di Simulazione Numerica  #
#                                       #
#             Lezione 05                #
#                                       #
#          Francesco Bollati            #
#                                       #
########################################*/

#include<iostream>
#include<fstream>
#include<string>
#include<cmath>

#include"random.h"
#include"funzioni.h"

using namespace std;


int main( int argc, char* argv[] ){

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



   //montecarlo simulation variables
   int teq;                                //equilibration time
   int Nsteps = 100000;                    //number Markov chain time steps 
   int step_block = 200;                   //number steps in a block 
   int N_blocks = int(Nsteps/step_block);  //number of blocks

   //DataBlock variables
   double sample[Nsteps];        
   double mean[N_blocks];
   double MeanError[N_blocks];


   //Markov chain starting points:
   //at argmax(|psi_100|^2)
   punto R0; 
   R0.x = 0;
   R0.y = 0;
   R0.z = 0;
   //far from argmax(|psi_100|^2)
   punto R_far;
   R_far.x = 10;
   R_far.y = 10;
   R_far.z = 10;
 
   punto R_lobo;
   R_lobo.x = 0;
   R_lobo.y = 0;
   R_lobo.z = 2;

   //Markov chain step length
   double r100; 
   double r210;


///////////////////////////////////////////// 5 0 %   A C C E P T A N C E /////////////////////////////////////////////

  
   //mettere commento in cui spiego cosa faccio

   punto Precedente;
   punto Attuale = R0;
   double m = 0;         //number accepted moves  

   for(int t = 0; t<Nsteps; t++){

      Precedente = Attuale; 
      Metropolis1(Attuale, 1, rnd, psi2_100); //generate next point 
      if( !(Precedente==Attuale) ) m++;       //accepted
      
    }
    //cout<<m/Nsteps<<endl;
    
   r100 = 1; 
   r210 = 2.5; 

///////////////////////////////////////// C O R R E L A T I O N   L E N G T H /////////////////////////////////////////


   
   ofstream out;
   out.open("correlazione.dat");

   int Omega = 10000;                //number of Markov chains
   punto R1[Omega];                  //array of points of time index = 1 for each Markov chain
   int tau_fin = 150;                //final time

   double Cov[tau_fin-1];            //cov(t=1,t=tau>1) with tau from tau = 2 to tau = taufin ( Cov[i] = cov(t=1,t=i+2) ) 
 
   for(int t =0; t<tau_fin -1; t++) Cov[t] = 0;

   for(int w = 0; w < Omega; w++){   //cycle on Markov chains

        //generete point R(t=1) for Markov chain w
        R1[w] = R0;
        Metropolis1(R1[w], r210, rnd, psi2_210);
        punto P = R1[w];

   	for(int t = 0; t < tau_fin -1 ; t++){

        	Metropolis1(P,r210,rnd, psi2_210);       
                Cov[t] += R1[w].dot(P)/double(Omega); 
                
                if(w == Omega -1) out<<t+2<<"   "<<Cov[t]<<endl; 
        }        
    }

    out.close();

    //fitting cov(1,tau) with a*exp(-x/l)+b I found l=8.4 for psi100 and l=34.4 for psi210. Then I set corrlation length at
    //least five times l. I choose step_block=200 for both cases. 



///////////////////////////////////////// E Q U I L I B R A T I O N   T I M E /////////////////////////////////////////



   out.open("equilibrio.dat");
   punto R;
 
   //R = R0;
   //R = R_far; 
   R = R_lobo;

   for(int t=0; t<500; t++){                     
   	Metropolis1(R,r210,rnd,psi2_210);
        out<<R.r()<<"   "<<R.z<<endl;
   }
   out.close();


////////////////////////////////////// M O N T E C A R L O   S I M U L A T I O N //////////////////////////////////////


   teq = 50; //starting from R0
   char* file = "psi2_210.dat"; //la string è const come in python quindi bisogna salvarlo in un CONST char*
   MonteCarlo(Nsteps, teq, Metropolis1, psi2_210, N_blocks, R0, rnd, file, r210, sample, mean, MeanError);


   teq = 50; //starting from R0
   char* file2 = "psi2_210_gauss.dat";
   //MonteCarlo(Nsteps, teq, Metropolis1, psi2_210, N_blocks, R0, rnd, file2, r210, sample, mean, MeanError);


   return 0;
}

