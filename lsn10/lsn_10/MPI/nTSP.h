#ifndef __GA_h__
#define __GA_h__

using namespace std;


#include "random.h"
int seed[4];
Random rnd; 

//simulation variables
int Ngenerations, Nchromosomes, Ncities, Noffspring, Ncores;
bool b, ll;
double pc;
double eps = 0.01;

//MPI
int rankm,size, indbest;
double bestL[100], lmin;


//genes, chromosomes, population
struct chromosome{
       
        double fitness;
        double l;

        int n_genes;
        int *genes;
        void equal( chromosome X ){
        
           l = X.l;
           n_genes = X.n_genes;
           fitness = X.fitness;
           genes = new int[n_genes];
           for(int i=0; i<n_genes; i++) genes[i] = X.genes[i];

        }
        void inizialize(int n){

          fitness = 0;
          l = 0;
          n_genes = n;
          genes = new int[n];
          

        }
};

chromosome population[1000]; 
chromosome NewPopulation[1000];

//temperature
double beta;
double beta_max,beta_min;
int n_beta;
double betav[1000]; 
int n_gen_per_T;

double accepted,attempted;


//city
struct city{
   
    double x;
    double y;

};

city *cities;


//functions
void Input(void);
void CreatePopulation(void);
void CalculateFitness(void);;
void Crossover(chromosome,chromosome,chromosome,chromosome);
void Mutate(void);
void ReplacePopulation(chromosome[]);
void Swap( chromosome, double );
void Check( chromosome );
void GenerateCities(bool);
void CalculateFitenss(void);
void QuickSort(chromosome[], int, int);
chromosome Selection(int&);
void Shift(chromosome,int);
void Swap_m( chromosome );
void Inversion( chromosome);
void Mutation( chromosome[]); 
void Metropolis(chromosome, chromosome, chromosome, chromosome);
void ReadCities(bool);
void Deallocate(void);
void PrintFinalSequence(void);
int SortPrint(int);
void CalculateL(double&, chromosome);

double L1( city[] );
double L2( city[] );





#endif
