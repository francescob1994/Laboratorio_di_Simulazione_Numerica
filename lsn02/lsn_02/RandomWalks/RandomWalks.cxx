/*#######################################
#                                       #
#  Laboratorio di Simulazione Numerica  #
#                                       #
#             Lezione 02                #
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

struct punto{ double x,y,z; };
punto incremento_reticolo(punto, Random&, double);   //generates i+1-th point given the i-th point in a 3D lattice
punto incremento_continuo(punto, Random&, double);   //generates i+1-th point given the i-th point in the continuum
void GenerateTrajectories( punto**, int, int, punto (*incremento)(punto, Random&, double), Random, double); //generates sample of trajectories 
void PrintTraj0( punto**, char*, int);

int main (int argc, char *argv[]){

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



   int M = 100;                  //Number of time step
   int n_traj = 10000;           //Number of trajectories
   double a = 1;                 //space increment


   //RANDOM WALK ON 3D LATTICE
   
   punto **R;                    //pointer to a pointer to allocate matrix for trajectories
   R = new punto *[n_traj];       
   for(int w=0; w <n_traj; w++){
   	R[w] = new punto[M];	
   }
  
   for(int w=0; w<n_traj; w++){  //set R = (0,0,0) in t = 0 for every trajectory
   	R[w][0].x = 0;
   	R[w][0].y = 0; 
   	R[w][0].z = 0;
   }
 
   GenerateTrajectories(R, M, n_traj, incremento_reticolo, rnd, a);

   char * file = "lattice.dat";
   PrintTraj0(R,file, M);




   //RANDOM WALK IN THE CONTINUUM

   punto **Rc;                    
   Rc = new punto *[n_traj];      
   for(int w=0; w <n_traj; w++){
   	Rc[w] = new punto[M];	
   }
  
   for(int w=0; w<n_traj; w++){  
   	Rc[w][0].x = 0;
   	Rc[w][0].y = 0; 
   	Rc[w][0].z = 0;
   }
 
   GenerateTrajectories(Rc, M, n_traj, incremento_continuo, rnd, a);

   file = "continuum.dat";
   PrintTraj0(Rc,file,M);




   //VARIANCE

   double Var[M-1];            //Variance_of_R = < |Rt(w)|^2 > for every t ( < Rt(w) >^2 = 0 ), lattice RW
   double VarError[M-1];       

   double Var_c[M-1];          //Variance_of_Rc, continuum RW
   double VarError_c[M-1];

   for(int t=0; t<M-1; t++){

      Var[t] = 0;
      double tmp = 0;
      double R2;

      Var_c[t] = 0;
      double tmp_c = 0;
      double Rc2;

      for(int w=0; w<n_traj; w++){

         R2 = R[w][t+1].x*R[w][t+1].x + R[w][t+1].y*R[w][t+1].y + R[w][t+1].z*R[w][t+1].z;
         Var[t] += R2;
         tmp += R2*R2;

         Rc2 = Rc[w][t+1].x*Rc[w][t+1].x + Rc[w][t+1].y*Rc[w][t+1].y + Rc[w][t+1].z*Rc[w][t+1].z;
	 Var_c[t] += Rc2;
         tmp_c += Rc2*Rc2; 
      }

      Var[t] /= n_traj;
      Var_c[t] /= n_traj;

      VarError[t] = sqrt( 1./(n_traj - 1) * ( tmp/n_traj - Var[t]*Var[t] ) );
      VarError_c[t] = sqrt( 1./(n_traj - 1) * ( tmp_c/n_traj - Var_c[t]*Var_c[t] ) );
   }

   
   //print data on diffiusion.dat
   
   ofstream diff;
   diff.open("diffusion.dat");
   for(int i=0; i<M-1; i++){
      diff<<i+1<<"   "<<Var[i]<<"   "<<VarError[i]<<"   "<<Var_c[i]<<"   "<<VarError_c[i]<<endl;
   }
   diff.close();   


   //deallocate
   for(int w=0; w<n_traj; w++){
      delete[] R[w];
      delete[] Rc[w];
   }




   rnd.SaveSeed();
   return 0;
}


//FUNCTIONS


punto incremento_reticolo( punto P, Random& r, double a ){ 
        punto R = P;   
        double y = r.Rannyu();
	if(y<1./6){
        	R.x+=a;	
        }else if(y>= 1./6 && y<1./3){
        	R.x-=a;        
        }else if(y>= 1./3 && y<1./2){
                R.y+=a;
        }else if(y>= 1./2 && y<2./3){
                R.y-=a;
        }else if(y>= 2./3 && y<5./6){
                R.z+=a;
        }else if(y>= 5./6){
                R.z-=a;
        }
        return R;
};

punto incremento_continuo( punto P, Random& r ,double a){
	punto R = P; 
        //genero phi
        double x = 1;
        double y = 1; 
        while(x*x+y*y>1){
        	x = r.Rannyu(-1,1);
                y = r.Rannyu(-1,1);
        }
        //genero theta
        double rho = 1;
        double z = 1;
        while(z*z+rho*rho >1){
        	rho = r.Rannyu();
                z = r.Rannyu(-1,1);
         }
         R.x += rho/sqrt(rho*rho + z*z) * x/sqrt(x*x + y*y);
         R.y += rho/sqrt(rho*rho + z*z) * y/sqrt(x*x + y*y);
         R.z += z/sqrt(rho*rho + z*z);
         return R;
};

void GenerateTrajectories( punto **R, int Nstep, int Ntraj, punto (*incremento)( punto P, Random& r, double a), Random r, double a){
	for(int w=0; w<Ntraj; w++){
		for(int t=0; t<Nstep-1; t++){
                	R[w][t+1] = (*incremento)(R[w][t],r,a);
		}
	}
}; 

void PrintTraj0( punto **R, char* file, int M){

   ofstream w;
   w.open(file);
   for(int i=0; i<M; i++){
      
      w << R[0][i].x << "   " << R[0][i].y << "   " << R[0][i].z << endl;

   }
   w.close();

}


