#ifndef __GA_h__
#define __GA_h__

using namespace std;


//simulation variables
int Ngenerations, Nchromosomes, Ncities, Noffspring;
bool b,ll;
double pc;
double eps = 0.01;


//genes, chromosomes, population
struct chromosome{
       
        double fitness;
        double l;      //length of correspondig path

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
void Swap(chromosome, double);
void Check(chromosome);
void GenerateCities(bool);
void ReadCities(bool);
void CalculateFitenss(void);
void QuickSort(chromosome[], int, int);
chromosome Selection(int&);
void Shift(chromosome,int);
void Swap_m(chromosome);  //swaps two blocks of m genes each
void Inversion(chromosome);
void Mutation(chromosome[]); 
void Deallocate(void);
void PrintLength(int);
void PrintFinalSequence(void);

double L1( city[] );
double L2( city[] );



#endif
