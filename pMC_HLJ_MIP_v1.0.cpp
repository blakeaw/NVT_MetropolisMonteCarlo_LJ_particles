/* NVT Metropolis Monte Carlo simulation of homogenous Lennard-Jones (LJ) particle system
Configured to run for reduced lj units with sigma=1.0 and epsilon=1.0 
using the standard 12-6 LJ potential for pair-wise interactions.
The Monte Carlo trial move is a single atom translation in 3d 
 The boundary condition is a rectangular box with min image periodicity; 
---Contains conditions to dynamically modulate the trial move step size.

Modifications


Author:Blake Wilson  
email:blake.wilson@utdallas.edu
Please e-mail me if you have comments, suggestions, or to report errors/bugs.
-------------------------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <time.h>
#include <sstream>

//lennard-jones cutoff distance
#define Rcut 2.5

using namespace std;

// Classes -- compiles fast enough to not worry about separately compiling classes
// Note-Since the classes are relatively small the .h file contains all the functions. 
// No additional .cpp files are used.
#include "class/Atom.h"
#include "class/Config.h"
#include "class/RunningStats.h"
#include "class/MTRNG.h"



// prototype functions

long double PairEnergy(const Atom& one, const Atom& two, const long double& boxx, const long double& boxy, const long double& boxz);
long double PairEnergyO(const Atom& one, const long double& ox, const long double& oy, const long double& oz, const long double& boxx, const long double& boxy, const long double& boxz);
long double TotalEnergy(const Config& c);
long double DeltaE(const Config& c, const unsigned& Aindex, const long double& ox, const long double& oy, const long double& oz);
string itos( int Number );
long double RDist(const Atom& one, const Atom& two);
long double RDistMinImage(const Atom& one, const Atom& two, long double bx, long double by, long double bz);
long double ConfProb(long double E, long double Beta);
string ftos (float input);

int main(int argc, char *argv[])
{

	//Get the local time for the logfile
	time_t rawtime;
 	struct tm * timeinfo;
 	time (&rawtime);
 	timeinfo = localtime (&rawtime);

	// Parameters that need to be set

	// the number of atoms to have in the system
	const int natoms = 25;
	
	//Total number of MC trial moves under each NS iteration
	const int MCit = 800000;

	// Number of initial MCit steps to devote to equilibration
	const int Eqit = 100000;

	// Step interval to add data to accumulated averages after equilibration
	const int Cit = natoms*2;
	
	//Step interval to output coordinates (.xyz) after equilibration
	unsigned int ConfigsOutFreq = ((MCit-Eqit))/20;

	//Temperature 
	long double Temp = 2.0;

	//Box sizes
	long double boxx = 5.0;
	long double boxy = 5.0;
	long double boxz = 5.0;
	
	//Initial step sizes
	long double stepx = boxx*0.5;
	long double stepy = boxy*0.5;
	long double stepz = boxz*0.5;

	//interval to check acceptance ratios and dynamically adjust
	//step sizes
	int ait = 4*natoms;
	
	//run descriptor
	string nds = itos(natoms);
	string ts = ftos(Temp);
	string rundescript = "MMC_n"+nds+"_c0_T"+ts+"_a";

	
//------------------------------------------------------------
	
	//Name of logfile
	string logfilename = "oOut_" + rundescript + ".log";
	//Coordinate output file
	string coutname = "cCoord_" + rundescript + ".xyz";
	
	ofstream logfile(logfilename.c_str());	
	ofstream coutfile(coutname.c_str());

	//check box sizes against lj cutoff
	long double trc = 2.0*Rcut;
	if (boxx<trc){
		cout<<"Warning! the x box size is less than two times the lj cutoff distance!"<<endl;
		cout<<"-Periodic interactions will not behave correctly. Adjust box size or cutoff."<<endl;
		exit(0); 
	}
	else if (boxy<trc){
		cout<<"Warning! the y box size is less than two times the lj cutoff distance!"<<endl;
		cout<<"-Periodic interactions will not behave correctly. Adjust box size or cutoff."<<endl;
		exit(0); 
	}
	if (boxz<trc){
		cout<<"Warning! the z box size is less than two times the lj cutoff distance!"<<endl;
		cout<<"-Periodic interactions will not behave correctly. Adjust box size or cutoff."<<endl;
		exit(0); 
	}

cout <<" MC iterations set to " << MCit << endl;
cout<<"Running at Temperature: "<<Temp<<" with "<<natoms<<" lj particles "<<endl;
logfile << " " << endl;
logfile << "Local time and date: " << asctime(timeinfo) << endl;
logfile<<"NVT MMC simulation of homogenous lennard jones system"<<endl;
logfile << "Parameters: " << endl;
logfile<<"Temperature: "<<Temp<<endl;
logfile<<" Number of Trial Moves (MCit): "<<MCit<<endl;
logfile<<"natoms: "<<natoms<<" Eqit: "<<Eqit<<" Cit: "<<Cit<<endl;
logfile<<" boxx: "<<boxx<<" boxy: "<<boxy<<" boxz: "<<boxz<<endl;
logfile<<"initial guess stepx: "<<stepx<<" stepy: "<<stepy<<" "<<stepz<<endl;
logfile<<"ConfigsOutFreq: "<<ConfigsOutFreq<<endl;

// Initialize the MT RNG object
MTRandomNum rng;


//Define a seed- use time
unsigned long long seed = time(0);
logfile<<"RNG seed: "<<seed<<endl;
logfile<<endl;
// Initialize the RNG object with the seed

rng.initialize(seed);
//Initialize running stat objects
//Potential energy
RunningStats Ustat;
//Potential energy squared -- for heat capacity
RunningStats Usqstat;

	

//Initialize configuration objects
	
Config trial1(natoms);
trial1.SetBox(boxx, boxy, boxz);


cout<<"Assigning Initial coordinates.."<<endl;
logfile<<"Assigning Initial coordinates.."<<endl;
//Randomly place atoms in box 

long double Et1;
long double Beta = 1.0/Temp;

//Assign a random initial configuration
for (int i = 0; i < natoms; ++i)
{
	trial1.atom[i].x = (rng.generate()-0.5)*boxx;
	trial1.atom[i].y = (rng.generate()-0.5)*boxy;
	trial1.atom[i].z = (rng.generate()-0.5)*boxz;
}
Et1 = TotalEnergy(trial1);
cout<<"Initial coordinates assigned with potential energy: "<<Et1<<endl;
logfile<<"Initial coordinates assigned with potential energy: "<<Et1<<endl;



//trial move acceptance trackers
int ctries = 0;
int csuccess = 0;
	
	
//Perform markov moves on the configuration
for(int n = 0; n<MCit; ++n){
		
		

	//Trial Move -----
	// randomly select an atom index	
	int Arandom = (int)(rng.generate()*natoms);
	//translations
	long double nx, ny, nz;
	nx = trial1.atom[Arandom].x + (rng.generate()-0.5)*stepx;
	ny = trial1.atom[Arandom].y + (rng.generate()-0.5)*stepy;
	nz = trial1.atom[Arandom].z + (rng.generate()-0.5)*stepz;
	//wrap coordinates
	trial1.WrapCoord(nx, ny, nz);
	//compute the change in potential energy
	long double dE=DeltaE(trial1, Arandom, nx, ny, nz);


	++ctries;
		
	//Selection Criteria--Metropolis Criteria

	//compute the probability ration
	long double Prat = ConfProb(dE, Beta);
	//random probability
	long double Prand = rng.generate();

	if(Prand < Prat){
		//trial move success
		++csuccess;	
		//udpdate energy		
		Et1+=dE;
		//update to new coordinates
		trial1.atom[Arandom].x=nx;
		trial1.atom[Arandom].y=ny;
		trial1.atom[Arandom].z=nz;
	}
	
// check acceptance ratio and dynamically adjust if needed
// target 33-55% acceptance
	if ( (n%ait)==0 ){

		long double locrat =  (long double)(csuccess)/(long double)(ctries);

		if (locrat<0.2){
			stepx*=0.1;
			stepy*=0.1;
			stepz*=0.1;

		}
		else if (locrat<0.33){

			long double stepc = locrat/0.33;
			stepx*=stepc;
			stepy*=stepc;
			stepz*=stepc;


		}
		else if(locrat>0.55){

			long double stepc = locrat/0.55;
			stepx*=stepc;
			stepy*=stepc;
			stepz*=stepc;

		}
					
	}
	// accumulate averages
	if( (n>Eqit) && ( ((n-Eqit)%Cit)==0 ) ){

		Ustat.Push(Et1);
		Usqstat.Push(Et1*Et1);

	}
			
		


//Output coordinates of sample n
	if( (n>Eqit) && ( ((n-Eqit)%ConfigsOutFreq)==0 ) ) {
		
		coutfile << trial1.natoms << endl;
		coutfile <<"lj "<<trial1.natoms<<" U: "<<Et1<<endl;
		coutfile << setiosflags(ios::fixed) << setprecision(10);
		for(int hg = 0; hg<trial1.natoms; ++hg){
			coutfile<<"lj "<<trial1.atom[hg].x<<" "<<trial1.atom[hg].y<<" "<<trial1.atom[hg].z<<endl; 
		}
		


	}	
		

		
} //End Monte Carlo loop



		
long double Uavg = Ustat.Mean();
long double Usqavg = Usqstat.Mean();
long double Cv = (Usqavg - Uavg*Uavg)*(Beta/Temp);
cout<<" Average Energy is "<<Uavg<<" and Average Energy Squared: "<<Usqavg<<endl;
cout<<" Heat capacity is: "<<Cv<<endl;
logfile<<"Average Energy is "<<Uavg<<" and Average Energy Squared: "<<Usqavg<<endl;
logfile<<" Heat capacity is: "<<Cv<<endl;
//close output files
logfile.close();
coutfile.close();
//complete
logfile<<"Simulation Complete!"<<endl;
logfile<<"Output files are: "<<endl;
logfile<<"Configuration Samples: "<<coutname<<endl;
cout<<"Simulation Complete!"<<endl;
cout<<"Output files are: "<<endl;
cout<<"Logfile: "<<logfilename<<endl;
cout<<"Configuration Samples: "<<coutname<<endl;
cout<<endl;
cout<<endl;
cout << " Thank You and Have a Nice Day!" << endl;
cout<<endl;
cout << "... End of Line ..." << endl;
	return EXIT_SUCCESS;



}

// Define Functions

// Total potential energy calculator, uses only the lj pot
long double TotalEnergy(const Config& c){
	
	long double Etot = 0.0;
	unsigned int natoms = c.natoms;
	long double bx, by, bz;
	bx = c.boxx, by = c.boxy, bz = c.boxz;
		for(unsigned int pp = 0;pp<natoms-1;++pp){

			for(unsigned int uu = pp+1;uu<natoms;++uu){
				
				
				Etot += PairEnergy(c.atom[pp], c.atom[uu], bx, by, bz);
					
			}
		}
//	Etot *= 60.0;
	return(Etot);
}

// potential energy function for LJ pair
long double PairEnergy(const Atom& one, const Atom& two, const long double& boxx, const long double& boxy, const long double& boxz){
	
	
	long double dx, dy, dz, rs;
	long double potent, eps, sigma;
	eps = 1.0, sigma = 1.0;
	long double bx=boxx;
	long double by=boxy;
	long double bz=boxz;
//	eps = sqrt(one.eps*two.eps), sigma = (one.sig+two.sig)/2.0;
	long double rc =Rcut;
	long double rcs = rc*rc;
	dx = one.x - two.x;
	
	dy = one.y - two.y;
	
	dz = one.z - two.z;
	//Minimum image -- coordinates must be pre-wrapped
	if(abs(dx)> bx/2.0){
		dx = bx - abs(one.x) - abs(two.x);
	}
	if(abs(dy)>by/2.0){
		dy = by - abs(one.y) - abs(two.y);
	}
	if(abs(dz)>bz/2.0){
		dz = bz - abs(one.z) - abs(two.z);
	}
	
	rs= (pow(dx, 2))+(pow(dy, 2))+(pow(dz, 2));
	if (rs<rcs){
		potent = 4.0*eps*((pow(sigma/rs, 6))-(pow(sigma/rs, 3))) - 4.0*eps*((pow(sigma/rcs, 6))-(pow(sigma/rcs, 3)));
	}
	else{
		potent=0.0;
	}	
	
	return(potent);
}

//Uses the new position values (xo, yo, zo) from atom one 
long double PairEnergyO(const Atom& one, const long double& ox, const long double& oy, const long double& oz, const long double& boxx, const long double& boxy, const long double& boxz){
	
	
	long double dx, dy, dz, rs;
	long double potent, eps, sigma;
	long double bx=boxx;
	long double by=boxy;
	long double bz=boxz;
	eps = 1.0, sigma = 1.0;
	long double rc = Rcut;
	long double rcs = rc*rc;
	dx = one.x - ox;
	
	dy = one.y - oy;
	
	dz = one.z - oz;

	//Minimum image -- coordinates must be pre-wrapped 
	if(abs(dx)> bx/2.0){
		dx = bx - abs(one.x) - abs(ox);
	}
	
	if(abs(dy)>by/2.0){
		dy = by - abs(one.y) - abs(oy);
	}
	
	if(abs(dz)>bz/2.0){
		dz = bz - abs(one.z) - abs(oz);
	}

	rs= (pow(dx, 2))+(pow(dy, 2))+(pow(dz, 2));
	if (rs<rcs){
		potent = 4.0*eps*((pow(sigma/rs, 6))-(pow(sigma/rs, 3))) - 4.0*eps*((pow(sigma/rcs, 6))-(pow(sigma/rcs, 3)));
	}
	else{
		potent=0.0;
	}	
	return(potent);
}

// calculates the change in energy associated with one moved atom
long double DeltaE(const Config& c, const unsigned& Aindex, const long double& ox, const long double& oy, const long double& oz){
	long double Etot=0.0;	
	long double Etoto=0.0;
	long double dE;
	long double bx, by, bz;
	bx = c.boxx, by = c.boxy, bz = c.boxz;
		
	
	for (unsigned int i = 0; i < c.natoms; ++i)
	{
		if(i!=Aindex){
			Etot += PairEnergy(c.atom[Aindex], c.atom[i], bx, by, bz);
			Etoto += PairEnergyO(c.atom[i], ox, oy, oz, bx, by, bz);
					
		}	
	}
	dE = Etot - Etoto;
	//cout << "In DeltaE funct, dE is : " << dE << endl;
	return( -dE ); 
	

}

string itos( int Number ){
     ostringstream ss;
     ss << Number;
     return ss.str();
 }

//Calculates the distance between two atoms--
long double RDist(const Atom& one, const Atom& two){

	long double dx, dy, dz, r, rs;
	long double x1, y1, z1, x2, y2, z2;
	x1 = one.x, y1=one.y, z1=one.z;
	x2=two.x, y2=two.y, z2=two.z;
	
	dx = x1-x2;
	
	dy = y1-y2;
	
	dz = z1-z2;
	
	rs= (pow(dx, 2))+(pow(dy, 2))+(pow(dz, 2));
	r = sqrt(rs);
	//cout<<"dx: "<<dx<<" dy: "<<dy<<" dz: "<<dz<<" r "<<r<<endl;
	return(r);



}
long double RDistMinImage(const Atom& one, const Atom& two, long double bx, long double by, long double bz){

	long double dx, dy, dz, r, rs;
	long double x1, y1, z1, x2, y2, z2;
	x1 = one.x, y1=one.y, z1=one.z;
	x2=two.x, y2=two.y, z2=two.z;
	
	dx = x1-x2;
	
	dy = y1-y2;
	
	dz = z1-z2;
	//Minimum image -- coordinates must be pre-wrapped
	if(abs(dx)> bx/2.0){
		dx = bx - abs(one.x) - abs(two.x);
	}
	if(abs(dy)>by/2.0){
		dy = by - abs(one.y) - abs(two.y);
	}
	if(abs(dz)>bz/2.0){
		dz = bz - abs(one.z) - abs(two.z);
	}
	rs= (pow(dx, 2))+(pow(dy, 2))+(pow(dz, 2));
	r = sqrt(rs);
	//cout<<"dx: "<<dx<<" dy: "<<dy<<" dz: "<<dz<<" r "<<r<<endl;
	return(r);



}


long double ConfProb(long double E, long double Beta){

	 
	return(exp(-Beta*E));
}

string ftos (float input){
	ostringstream strs;
	strs << input;
	string str = strs.str();
	return str;

}
