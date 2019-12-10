/*#######################################
#                                       #
#  Leboratorio di Simulazione Numerica  #
#                                       #
#             Lezione 10                #
#                                       #
#          Francesco Bollati            #
#                                       #
########################################*/


#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<stdlib.h>
#include"TSP.h"
#include <time.h>


using namespace std;

int main()
{

   Input();

   clock_t start,end;
   double tempo;
   start=clock();
   
   CreatePopulation();
   CalculateFitness();
     
   for(int i=0; i<n_beta; i++){  //cycle on temperatures

   beta = betav[i]; 

   attempted = 0;
   accepted = 0;

      for( int t = 1; t <= n_gen_per_T; ++t) //cycle on generation per temperature
      {

         Noffspring = 0;

         while( Noffspring < Nchromosomes )
         {
            //SELECTION OF PARENTS

            int i,j;                      //indices of selected chromosomes
            chromosome X = Selection(i); 
            chromosome Y = Selection(j);

            int exit=0;
            while( i == j ){
               if(exit==100){
                  cerr<<"It can't select a parent Y different from X, too small probability! Try again!"<<endl;
                  return -1;
               }
               delete[] Y.genes;
               Y.equal(Selection(j)); 
            }

            //CROSSOVER, ALLOCATE MEMORY FOR OFFSPROING

            if(drand48()<pc)
            {

               chromosome XX,YY;         //define offsprings and allocate memory for them
               XX.inizialize(X.n_genes);
               YY.inizialize(X.n_genes);

               Crossover(X,Y,XX,YY);    //fills offsprings' genes

               Metropolis(X,Y,XX,YY);   //updates NewPopulation and Noffspring

               delete[] XX.genes; 
               delete[] YY.genes; 

            }else{  

               NewPopulation[Noffspring].equal(X);     
               NewPopulation[Noffspring+1].equal(Y);
               Noffspring+=2;
               delete[] X.genes; 
               delete[] Y.genes;  

            }
               
         }

         Mutation(NewPopulation); 

         ReplacePopulation(NewPopulation); 

         CalculateFitness();

         PrintLength(i*n_gen_per_T + t);

      }

      //cout<<"beta = "<<beta<<": accepted/attempted = "<<accepted/attempted<<endl;
   }

   PrintFinalSequence();

   Deallocate(); 

   end=clock();
   tempo=((double)(end-start))/CLOCKS_PER_SEC;
   cout<<"tempo di esecuzione "<<tempo<<endl;  

return 0;
}


void Input(void)
{
   cout << "Solution of the Traveling Salesman Probelm with Genetic Algorithm and Simulated Annealing" << endl <<endl;

   srand48(time(NULL));
  
   ifstream ReadInput;
   ReadInput.open("input.dat");
   
   ReadInput >> Ncities;
   ReadInput >> Nchromosomes;
   ReadInput >> n_gen_per_T;
   ReadInput >> b;
   ReadInput >> pc;
   ReadInput >> ll;

   ReadInput >> n_beta;
   ReadInput >> beta_max;
   ReadInput >> beta_min;

   cout << "Number of cities = " << Ncities << endl;
   cout << "Number of chromosomes = " << Nchromosomes << endl;
   cout << "Number of generations per fixed Temperature = " << n_gen_per_T <<endl;
   cout << "Number of different temperatures = "<< n_beta <<endl;
   cout << "Beta max = " << beta_max <<endl;
   cout << "Beta min = " << beta_min <<endl;

   //GenerateCities(b); //generate positions of cities
   ReadCities(b);

   //inizialize vector of betas;
   double delta_beta = (beta_max - beta_min)/double(n_beta - 1);
   for(int i = 0; i<n_beta; i++){
      betav[i] = beta_min + i*delta_beta; 
   }
 
}

void ReadCities(bool b)
{
   ifstream read;
   if(b) read.open("square.dat");
   else read.open("circumference.dat");

   cities = new city[Ncities];
   for(int i=0; i<Ncities; i++){
      read >> cities[i].x;
      read >> cities[i].y;
   }
   read.close();   

}


void CreatePopulation(void)
{
   chromosome X;
   X.genes = new int[Ncities];
   X.n_genes = Ncities;
   for(int i=0; i<Ncities; i++) X.genes[i] = i;

   for(int i=0; i<Nchromosomes; i++)
   {
      population[i].equal(X);
      Swap(population[i],0.2);
      Check( population[i] );
  
   }
}

void Check( chromosome X )
{
   int N = X.n_genes;
   for(int i=0; i<N-1; i++)
   {
      for(int j=i+1; j<N; j++)
      {
         if(X.genes[i] == X.genes[j])
         {
            cout<<"Repetition of genes"<<endl;
            exit(-1);
         }
      }
    }
}

void GenerateCities(bool b)
{

    cities = new city[Ncities];
    if(b){   //square
 
       double L = sqrt(Ncities);
       for(int i=0; i<Ncities; i++)
       {
          cities[i].x = 2*L*drand48() - L;
          cities[i].x = 2*L*drand48() - L;     
       }

    }else{   //circumference

       double R = sqrt(Ncities);
       for(int i=0; i<Ncities; i++)
       {
          double phi = drand48()*2*M_PI;
          cities[i].x = R*cos(phi);
          cities[i].x = R*sin(phi);
       }
    }
}

double L1( city cities[] )
{
   double dist = 0; 
   for(int i=0; i<Ncities-1; i++)
   {
      dist += sqrt(  pow(cities[i+1].x -cities[i].x,2) + pow(cities[i+1].y -cities[i].y,2)  );
   }
   dist += sqrt(  pow(cities[Ncities-1].x -cities[0].x,2) + pow(cities[Ncities-1].y -cities[0].y,2)  );
   return dist;
}

double L2( city cities[] )
{
   double dist = 0; 
   for(int i=0; i<Ncities-1; i++)
   {
      dist += pow(cities[i+1].x -cities[i].x,2) + pow(cities[i+1].y -cities[i].y,2);
   }
   dist += pow(cities[Ncities-1].x -cities[0].x,2) + pow(cities[Ncities-1].y -cities[0].y,2);
   return dist;
}


void CalculateFitness(void)
{  

   //calculate the length of chromosomes
   for(int i=0; i<Nchromosomes; i++)
   {
     CalculateL( population[i].l, population[i] );
   }

   //sort the population in order of descending fitness
   QuickSort(population, 0, Nchromosomes-1);


   //calculate the fitness of chromosomes
   double l_max = population[Nchromosomes-1].l;
   double sum = 0;
   for(int i = 0; i<Nchromosomes; i++)
   {
      population[i].fitness =-exp(-0.1*l_max)*(1-0.5) + exp(-0.1*population[i].l);//-population[i].l + l_max*(1+0.1); // 
      sum += population[i].fitness;        
 
   }

   for(int i = 0; i<Nchromosomes; i++)
   {
      population[i].fitness = population[i].fitness/sum;  //the normalized fitness is the probability to choose the chromosome
      if(population[i].fitness == 0 ) population[i].fitness = eps;
      if(population[i].fitness == 1 ) population[i].fitness = 1 - eps;

   }    
   
}


//sort in order of ascending length
void QuickSort( chromosome pop[], int low, int high ) 
{                                    
   int l = low;
   int r = high;
   double pivot = pop[(l+r)/2].l;    
   chromosome tmp;

   while(l<=r)
   {
      while(pop[l].l < pivot) l++;
      while(pop[r].l > pivot) r--;
      if( l<=r ){

         tmp.equal(pop[l]);
         delete[] pop[l].genes;
         pop[l].equal(pop[r]);
         delete[] pop[r].genes;
         pop[r].equal(tmp);
         delete[] tmp.genes;

         l++;
         r--;
      }
   }
   if(r>low) QuickSort( pop, low, r);
   if(l<high) QuickSort( pop, l, high); 
}


chromosome Selection(int& k)
{

   chromosome X;

   //create a partition of [0,1]: [0, population[0].fitness, population[1].fitness, .. ]
   double x[Nchromosomes+1];
   x[0] = 0;
   for(int i=1; i<Nchromosomes+1; i++)
   {
      x[i] = population[i-1].fitness + x[i-1];
   }


   double r = drand48();
   for(int i = 0; i<Nchromosomes; i++)
   {
      if( r < x[i+1] && r >= x[i] ){
         k = i;
         X.equal(population[i]);
         return X;       
         break;
      }
    }
}


//MUTATIONS

void Shift( chromosome X, int n )
{
   int N = X.n_genes;

   int v[N];

   for(int i=0; i<N; i++)
   {
      v[(n+i)%N] = X.genes[i];
   }
   for(int i=0; i<N; i++)
   {
      X.genes[i] = v[i];
   }
}

void Swap( chromosome X , double p) 
{
   int N = X.n_genes;
   int tmp;
   for(int i=0; i<N-1; i++)
   {
      for(int j=i+1; j<N; j++)
      {
         if(drand48()<p)
         {
            tmp = X.genes[i];
            X.genes[i] = X.genes[j];
            X.genes[j] = tmp;
         }
      }
    }
}

void Swap_m( chromosome X )
{   

        int N = X.n_genes;
        int m = 1+ int(drand48()*N/2); //random integer between 1 and N/2 (length of the blocks to swap)
   
        //select the first index of the left block such that starting from it there are at least 2m genes before the end of the chromosome
        int l = int(drand48()*(N-m*2+1));
        //select the first index of the right block such that it is grater than the last index of the first block and that starting
        //from it there are at least m genes before the end of the chromosome
        int r = int(l+m+ drand48()*(N-2*m-l));

        int v[m];
        for(int i=0; i<m; i++){
           v[i] = X.genes[l+i];
           X.genes[l+i] = X.genes[r+i];
           X.genes[r+i] = v[i];
        }

}


void Inversion( chromosome X )
{
        int N = X.n_genes;
        int m =  int(-(N/4)*log(1-drand48())); //integer drown from exponential distribution with mean = N/4

        if(m<=N && m>1){
        int l = int(drand48()*(N-m));
        double tmp;
        for(int i=0; i<m/2; i++)
        {    
           tmp = X.genes[l+i];
           X.genes[l+i] = X.genes[l+m-1-i];
           X.genes[l+m-1-i] = tmp;
        }
        }
}



void Mutation(chromosome NewPopulation[])
{
   int N = NewPopulation[0].n_genes;
   for(int i=0; i<Nchromosomes; i++)
   {

     int k = int(drand48()*5); //select one of the four possible mutations

     //once a mutation is selected it is applied with probability about 10%

     if(k == 0){

        if(drand48()<0.12){
          Shift(NewPopulation[i],2);
        }
   
     }else if( k == 1){

        double p = 0.12*2./(N*(N-1));  //probability of a single swap in the chromosome
        Swap(NewPopulation[i],p);
  
     }else if(k == 2){

        double p = 0.12/( exp(-6./double(N))-exp(-3.) ); //***************
        if(drand48()<p) Swap_m(NewPopulation[i]);
      
     }else{
     
        if(drand48()<0.12){
           Inversion(NewPopulation[i]); 

        }
     }
     
   }
}


void Crossover(chromosome X, chromosome Y, chromosome XX, chromosome YY)
{
   int N = X.n_genes;
   int ind = int( drand48()*N ); //index of slicing

   int n = 0, m = 0;                    
   int v[N-ind], w[N-ind];

   for( int i=0; i<N; i++)   
   {
      for( int j=ind; j<N; j++)  
      {
         if( X.genes[i] == Y.genes[j] )
         {
         v[n] = X.genes[i];
         n++;
         }
         if( Y.genes[i] == X.genes[j] )
         {
         w[m] = Y.genes[i];
         m++;
         }
      }
   }

           
 for( int i=0; i<N; i++)
   {
      if( i < ind )
      {
         XX.genes[i] = X.genes[i];
         YY.genes[i] = Y.genes[i];
      }
      else
      {
         XX.genes[i] = w[i-ind];
         YY.genes[i] = v[i-ind];
      }
    }

}


void ReplacePopulation( chromosome NewPopulation[] )
{
   for(int i = 0; i<Nchromosomes; i++)
   {
      delete[] population[i].genes;
      population[i].equal(NewPopulation[i]);  
      delete[] NewPopulation[i].genes;
   }
    
}

void CalculateL(double &l, chromosome X)
{
      int N = X.n_genes;
      city cit[N]; 
      for(int j=0; j<N; j++)  //cycle on i-th chromosome's genes
      {
         int ind = X.genes[j];
         cit[j] = cities[ind];
      }
      if(ll) l = L2(cit);
      else l = L1(cit);
   
}



void Metropolis(chromosome X, chromosome Y, chromosome XX, chromosome YY)
{

 //sort X and Y
 chromosome C1,C2;
 if(X.l < Y.l){
    C1.equal(X); 
    C2.equal(Y);
 }else{
    C1.equal(Y);
    C2.equal(X);
 }

 delete[] X.genes; 
 delete[] Y.genes; 

 bool used;                               //true if Metropolis rejects the first offspring and updates NewPopulation with parent C1

 CalculateL(XX.l, XX);                    //calculate length of the first offspring

 double p = exp(-beta*XX.l)/exp(-beta*C1.l); 

 if(p>=drand48()){

    NewPopulation[Noffspring].equal(XX);  //accept the first offspring
    accepted ++; 
    used = false;

  }else{

    NewPopulation[Noffspring].equal(C1);  //reject the first offspring and update NewPopulation with the parent with the best fitness
    used = true;

  } 

  Noffspring++;  
  attempted++;

 CalculateL(YY.l, YY);                    //calculate length of the second offspring
 p = exp(-beta*YY.l)/exp(-beta*C1.l); 

 if(p>=drand48()){

    NewPopulation[Noffspring].equal(YY);  //accept the second offspring
   accepted++;

 }else{
   
    if(!used){                            //if the second offspring is rejected and C1 hasn't been used in the previous update
                                          //C1 is added to NewPopulation
       NewPopulation[Noffspring].equal(C1);

    }else{                                //if the second offspring is rejected and C1 has been used in the previous update
                                          //C2 is added to NewPopulation
       NewPopulation[Noffspring].equal(C2);
  
    }
  } 

  Noffspring++;
  delete[] C1.genes; 
  delete[] C2.genes;
  attempted ++;
}



void Deallocate(void)
{
   for(int i = 0; i<Nchromosomes; i++)
   {
      delete[] population[i].genes;
   }

}

void PrintLength(int i)
{
   ofstream len;
   len.open("length.dat", ios::app);

   double l_ave=0;
   for(int k=0; k<Nchromosomes/2; k++) l_ave += population[k].l;
   l_ave /= Nchromosomes/2;
   
   len << i << "   " << population[0].l << "   " << l_ave << endl;

   len.close();

}

void PrintFinalSequence(void)
{
   ofstream opt;
   opt.open("FinalSeq.dat");
   opt << b <<endl;
   opt << ll <<endl;
   opt << population[0].l << endl;
   for(int i=0; i<Ncities; i++){
      opt << population[0].genes[i] <<endl;
   }
   opt.close();

}




















