/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"
#include <iomanip>

using namespace std;

int main(){ 
  Input();             //Inizialization
  int nconf = 1;
  
  for(int iblk=1; iblk<=nblocks; iblk++){

     Reset(iblk);

     for(int istep=1; istep <= StepsPerBlock; ++istep){

        Move();        //Move particles with Verlet algorithm
        if( ((iblk-1)*StepsPerBlock+istep)%iprint == 0 ) cout << "Number of time-steps: " << ((iblk-1)*StepsPerBlock+istep) << endl;
        if(istep%10 == 0){
           Measure();  //Properties measurement
           //ConfXYZ(nconf);
           Accumulate();
           nconf += 1;          
         }
     }
     MeasureAve(iblk);
  }


  ConfOld();           //Write second-last configuration
  ConfFinal();         //Write final configuration to restart

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> nblocks;
  ReadInput >> old;

  StepsPerBlock = nstep/nblocks;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl;
  cout << "Number of blocks = " << nblocks << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//Read initial configuration

  if(old){

     //Read from old.final
     cout << "Read old configuration from file old.final " << endl << endl;
     ReadConf.open("old.final");
     for (int i=0; i<npart; ++i){
        ReadConf >> xold[i] >> yold[i] >> zold[i];
        xold[i] = xold[i] * box;
        yold[i] = yold[i] * box;
        zold[i] = zold[i] * box;
     }
     //Read from config.final
     ReadConf.close();
     cout << "Read initial configuration from file config.final " << endl << endl;
     ReadConf.open("config.final");
     for (int i=0; i<npart; ++i){
        ReadConf >> x[i] >> y[i] >> z[i];
        x[i] = x[i] * box;
        y[i] = y[i] * box;
        z[i] = z[i] * box;
     }
     ReadConf.close();

     //Create new configuration: xold, x ---> xnew ---> xold=x, x=xnew
     Move(); 
    
     //Find v(t+dt/2) with x, xold
     double sumv2 = 0.0, fs;
     for (int i=0; i<npart; ++i){
        vx[i] = Pbc(x[i]-xold[i])/delta; 
        vy[i] = Pbc(y[i]-yold[i])/delta;
        vz[i] = Pbc(z[i]-zold[i])/delta;

        sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
     }
     sumv2 /= (double)npart;

     //riscale velocities to get target temperature
     fs = sqrt(3 * temp / sumv2);    
     for (int i=0; i<npart; ++i){
        vx[i] *= fs;
        vy[i] *= fs;
        vz[i] *= fs;

        xold[i] = Pbc(x[i] - vx[i] * delta);
        yold[i] = Pbc(y[i] - vy[i] * delta);
        zold[i] = Pbc(z[i] - vz[i] * delta);
     }
     



  }else{

     cout << "Read initial configuration from file config.0 " << endl << endl;
     ReadConf.open("config.0");
     for (int i=0; i<npart; ++i){
        ReadConf >> x[i] >> y[i] >> z[i];
        x[i] = x[i] * box;
        y[i] = y[i] * box;
        z[i] = z[i] * box;
     }
     ReadConf.close();

     //Prepare initial velocities
     cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
     double sumv[3] = {0.0, 0.0, 0.0};
     for (int i=0; i<npart; ++i){
        vx[i] = rand()/double(RAND_MAX) - 0.5;
        vy[i] = rand()/double(RAND_MAX) - 0.5;
        vz[i] = rand()/double(RAND_MAX) - 0.5;

        sumv[0] += vx[i];
        sumv[1] += vy[i];
        sumv[2] += vz[i];
     }
     for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
     double sumv2 = 0.0, fs;
     for (int i=0; i<npart; ++i){
        vx[i] = vx[i] - sumv[0];
        vy[i] = vy[i] - sumv[1];
        vz[i] = vz[i] - sumv[2];

        sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
     }
     sumv2 /= (double)npart;

     fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
     for (int i=0; i<npart; ++i){
        vx[i] *= fs;
        vy[i] *= fs;
        vz[i] *= fs;

        xold[i] = Pbc(x[i] - vx[i] * delta);
        yold[i] = Pbc(y[i] - vy[i] * delta);
        zold[i] = Pbc(z[i] - vz[i] * delta);
     }

  }


  //g(r)
  nbins = 100;
  bin_size = (box/2.0)/(double)nbins;
  for(int i=0; i<nbins; i++) walker[i]=0.0;
  

  return;

}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij, w, wij; //
  double dx, dy, dz, dr;

  ofstream Epot, Ekin, Etot, Temp, Pres; //

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);
  Pres.open("output_pres.dat",ios::app); 

  v = 0.0; //reset observables
  t = 0.0;
  w = 0.0;

  for (int k=0; k<nbins; ++k) walker[k]=0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

//g(r)

     for(int k=0; k<nbins; k++){
        if(dr >= k*bin_size && dr < (k+1)*bin_size ) walker[k]++;
     }

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       wij = 48.0/pow(dr,12) - 24.0/pow(dr,6); 
//Potential energy
       v += vij;
       w += wij; 
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle
    stima_pres = rho*stima_temp + 1./(3.*vol) * w; //

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Pres << stima_pres << endl; //   

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Pres.close(); //

}

void Accumulate(){

  v_blk += stima_pot;//
  //v_ave2 += stima_pot*stima_pot;//

  k_blk += stima_kin;//
  //k_ave2 += stima_kin*stima_kin;//

  t_blk += stima_temp;//
  //t_ave2 += stima_temp*stima_temp;//

  e_blk += stima_etot;//
  //e_ave2 += stima_etot*stima_etot;//

  p_blk += stima_pres; //
  //p_ave2 += stima_pres*stima_pres; //

  n_mis++;
  //cout<<n_mis<<"  "<<k_blk<<endl;
  
  //g(r)
  for(int i=0; i<nbins; i++){
    blk_gofr[i] += walker[i];
  }
  
}


void MeasureAve(int iblk){

  ofstream EpotAve, EkinAve, EtotAve, TempAve, PresAve;

  EpotAve.open("ave_epot.out",ios::app);
  EkinAve.open("ave_ekin.out",ios::app);
  TempAve.open("ave_temp.out",ios::app);
  EtotAve.open("ave_etot.out",ios::app);
  PresAve.open("ave_pres.out",ios::app);
  
  v_blk = v_blk/double(n_mis);
  v_prog += v_blk;
  v_prog2 += v_blk*v_blk;
  v_err = Error(v_prog,v_prog2,iblk);

  k_blk = k_blk/double(n_mis);
  k_prog += k_blk;
  k_prog2 += k_blk*k_blk;
  k_err = Error(k_prog,k_prog2,iblk);

  //cout<<iblk<<" "<<k_blk<<"  "<<k_prog/double(iblk)<<"  "<<k_err<<endl;

  t_blk = t_blk/double(n_mis);
  t_prog += t_blk;
  t_prog2 += t_blk*t_blk;
  t_err = Error(t_prog,t_prog2,iblk);

  e_blk = e_blk/double(n_mis);
  e_prog += e_blk;
  e_prog2 += e_blk*e_blk;
  e_err = Error(e_prog,e_prog2,iblk);

  p_blk = p_blk/double(n_mis);
  p_prog += p_blk;
  p_prog2 += p_blk*p_blk;
  p_err = Error(p_prog,p_prog2,iblk);


  EpotAve << iblk << "    " << v_blk << "    " << v_prog/(double)iblk << "    " << v_err << endl;
  EkinAve << iblk << "    " << k_blk << "    " << k_prog/(double)iblk << "    " << k_err << endl;
  TempAve << iblk << "    " << t_blk << "    " << t_prog/(double)iblk << "    " << t_err << endl;
  EtotAve << iblk << "    " << e_blk << "    " << e_prog/(double)iblk << "    " << e_err << endl;
  PresAve << iblk << "    " << p_blk << "    " << p_prog/(double)iblk << "    " << p_err << endl;
	

  EpotAve.close();
  EkinAve.close();
  TempAve.close();
  EtotAve.close();
  PresAve.close();



   //g(r)
   ofstream Gofr,Gave;
   Gofr.open("gofr.dat",ios::app);
   Gave.open("gofr_ave.dat",ios::app);
   int wd=12;

   Gofr << iblk;
   double er;
   double dV;
   double norma;
   for(int k = 0; k<nbins; k++){
      er = k*bin_size;
      dV = (4*pi/3.)*( pow(er+bin_size,3)-pow(er,3) );
      norma = rho*npart*dV;

      stima_gofr = (blk_gofr[k]/(double)n_mis)/norma;
      glob_gofr[k] += stima_gofr;
      glob_gofr2[k] += stima_gofr*stima_gofr;
      err_gofr = Error(glob_gofr[k],glob_gofr2[k],iblk);

      Gofr << setw(wd) << er << setw(wd) << stima_gofr << setw(wd) << glob_gofr[k]/(double)iblk << setw(wd) << err_gofr << endl;

      if(iblk==nblocks){
         Gave << er << "   " << glob_gofr[k]/(double)nblocks << "   " << err_gofr <<endl; 
      }

   }



}


double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}



void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
   v_prog = 0;
   k_prog = 0;
   t_prog = 0;
   e_prog = 0;
   p_prog = 0;

   v_prog2 = 0;
   k_prog2 = 0;
   t_prog2 = 0;
   e_prog2 = 0;
   p_prog2 = 0;

   for(int i=0; i<nbins; i++) glob_gofr[i]=0;
   for(int i=0; i<nbins; i++) glob_gofr2[i]=0;
   
   }

   v_blk = 0;
   k_blk = 0;
   t_blk = 0;
   e_blk = 0;
   p_blk = 0;

   for(int i=0; i<nbins; i++) blk_gofr[i]=0;


   n_mis = 0;

}





void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}


void ConfOld(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print second-last configuration to file old.final " << endl << endl;
  WriteConf.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  return;
}


void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
