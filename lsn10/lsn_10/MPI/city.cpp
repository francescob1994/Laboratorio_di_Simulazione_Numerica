#include<iostream>
#include<fstream>
#include<cmath>
#include<stdlib.h>
#include<time.h>

using namespace std;

struct city{
   
    double x;
    double y;

};


int main(){

    srand48(time(NULL));

    ofstream ci;
    ci.open("città.dat");
   
    //quando valuto i tempi su numeri di città diverse devo aggiornare questo file e stmpare le nuove città
    int Ncities = 15; 
 
       //SQUARE
       double L = sqrt(Ncities);
       for(int i=0; i<Ncities; i++)
       {
           ci << (2*L*drand48() - L) <<"   "<<(2*L*drand48() - L)<< endl;

       }

      //CIRCLE
       double R = sqrt(Ncities);
     /*  for(int i=0; i<Ncities; i++)
       {
          double phi = drand48()*2*M_PI;
           ci << R*cos(phi) <<"   "<<R*sin(phi)<< endl;
       }*/  

       ci.close();
return 0;

}



