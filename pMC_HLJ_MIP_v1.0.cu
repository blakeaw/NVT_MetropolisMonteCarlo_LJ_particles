/* NVT Metropolis Monte Carlo simulation of homogenous Lennard-Jones (LJ) particle system
Configured to run for reduced lj units with sigma=1.0 and epsilon=1.0 
using the standard 12-6 LJ potential for pair-wise interactions.
The Monte Carlo trial move is a single atom translation in 3d 
 The boundary condition is a rectangular box with min image periodicity; 
---Contains conditions to dynamically modulate the trial move step size.

cuda parallelized for nvidia GPGPU

The inner loop of the total potential energy computation is parallized (only called
once in this code though)
and the potential energy difference computation O(2N-2) is parallized
-Requires an inclusive scan function: currently uses the Thrust 
(https://thrust.github.io/) library implementation for device
parallel inclusive_scan

Modifications


Author:Blake Wilson  
email:blake.wilson@utdallas.edu
Please e-mail me if you have comments, suggestions, or to report errors/bugs.

compile (nvidia Tesla compute 3.5):
nvcc -O2 -arch=compute_35 pMC_HLJ_MIP_v1.0.cu
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
//for cuda
#include <cuda.h>
//need the inclusive scan - use thrust
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
//lennard-jones cutoff distance
#define Rcut 2.5
//maximum thread per block for nvidia gpu
#define MAX_THREAD 256
using namespace std;

// Classes -- compiles fast enough to not worry about separately compiling classes
// Note-Since the classes are relatively small the .h file contains all the functions. 
// No additional .cpp files are used.
#include "class_cuda/Atom-b.h"
#include "class_cuda/Config-b.h"
#include "class_cuda/RunningStats.h"
#include "class_cuda/ConfigGPU-c.h"
#include "class_cuda/MTRNG.h"

// prototype functions
__device__ double PairEnergy_Cuda(double ax1, double ay1, double az1, double ax2, double ay2, double az2, double bx, double by, double bz);
__global__ void DEloop(unsigned int natoms, const unsigned int Aind, double * box, double * olc, double * ax, double * ay, double * az, double * dE);
__global__ void TEloop(int natoms, int Aindex, double *box, double *ax, double *ay, double *az, double *d_Ep, int sindex);
__global__  void GetDoubleDeviceArrayAt_Kernel(double *d_array, int index, double * d_v);
double GetDoubleDeviceArrayAt(double * d_array, int index);
__global__ void ZeroDeviceArray_Kernel(double *d_arr, int length);
void ZeroDeviceArray(double *d_array, int nelements);
double TotalEnergy_Cuda(const ConfigGPU& c);
double DeltaE_Cuda(const ConfigGPU& c, const unsigned& Aindex, const double& ox, const double& oy, const double& oz);
string itos( int Number );
double ConfProb(double E, double Beta);
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
	const int natoms = 62;
	
	//Total number of MC trial moves under each NS iteration
	const int MCit = 400000;

	// Number of initial MCit steps to devote to equilibration
	const int Eqit = 200000;

	// Step interval to add data to accumulated averages after equilibration
	const int Cit = 1000;
	
	//Step interval to output coordinates (.xyz) after equilibration
	unsigned int ConfigsOutFreq = ((MCit-Eqit))/2;

	//Temperature 
	double Temp = 2.0;

	//Box sizes
	double boxx = 5.0;
	double boxy = 5.0;
	double boxz = 5.0;
	
	//Initial step sizes
	double stepx = boxx*0.5;
	double stepy = boxy*0.5;
	double stepz = boxz*0.5;

	//interval to check acceptance ratios and dynamically adjust
	//step sizes
	int ait = 10000;
	
	//run descriptor
	string nds = itos(natoms);
	string ts = ftos(Temp);
	string rundescript = "MMC_n"+nds+"_c0_T"+ts+"_a";

	
//------------------------------------------------------------

//Availiable GPU data
int deviceCount;
cudaGetDeviceCount(&deviceCount);
if(deviceCount==0){
	cout<<"No nvidia GPUs found! exiting..."<<endl;
	exit(0);
}
int device;
int computemax=0;
int dcmax=0;
for (device = 0; device < deviceCount; ++device) {
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, device);
	printf("Device %d has compute capability %d.%d.\n",
	device, deviceProp.major, deviceProp.minor);
	int compute = 10*deviceProp.major + deviceProp.minor;
	if(compute>computemax){
		dcmax=device;
		computemax=compute;
	}
}
//set to one with highest compute capability
cout<<"setting device to :"<<dcmax<<" with compute "<<computemax<<endl;
cudaSetDevice(dcmax);



	//Name of logfile
	string logfilename = "oOut_" + rundescript + ".log";
	//Coordinate output file
	string coutname = "cCoord_" + rundescript + ".xyz";
	
	ofstream logfile(logfilename.c_str());	
	ofstream coutfile(coutname.c_str());

	//check box sizes against lj cutoff
	double trc = 2.0*Rcut;
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

	

//Initialize standard configuration object	
Config trial1(natoms);
trial1.SetBox(boxx, boxy, boxz);


cout<<"Assigning Initial coordinates.."<<endl;
logfile<<"Assigning Initial coordinates.."<<endl;
//Randomly place atoms in box 

double Et1;
double Beta = 1.0/Temp;

//Assign a random initial configuration
for (int i = 0; i < natoms; ++i)
{
	trial1.atom[i].x = (rng.generate()-0.5)*boxx;
	trial1.atom[i].y = (rng.generate()-0.5)*boxy;
	trial1.atom[i].z = (rng.generate()-0.5)*boxz;
}

//now initialize the GPU config object and copy the initial config over
ConfigGPU dconf(natoms);
dconf.FullEquateCPU(trial1);

Et1 = TotalEnergy_Cuda(dconf);

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
	double nx, ny, nz;
	nx = dconf.ax[Arandom] + (rng.generate()-0.5)*stepx;
	ny = dconf.ay[Arandom] + (rng.generate()-0.5)*stepy;
	nz = dconf.az[Arandom] + (rng.generate()-0.5)*stepz;
	//wrap coordinates
	dconf.WrapCoord(nx, ny, nz);
	//compute the change in potential energy
	double dE=DeltaE_Cuda(dconf, Arandom, nx, ny, nz);


	++ctries;
		
	//Selection Criteria--Metropolis Criteria

	//compute the probability ration
	double Prat = ConfProb(dE, Beta);
	//random probability
	double Prand = rng.generate();

	if(Prand < Prat){
		//trial move success
		++csuccess;	
		//udpdate energy		
		Et1+=dE;
		//update to new coordinates
		dconf.SetAtomCoord(Arandom, nx, ny, nz);
	}
	
// check acceptance ratio and dynamically adjust if needed
// target 33-55% acceptance
	if ( (n%ait)==0 ){

		double locrat =  (double)(csuccess)/(double)(ctries);

		if (locrat<0.2){
			stepx*=0.1;
			stepy*=0.1;
			stepz*=0.1;

		}
		else if (locrat<0.33){

			double stepc = locrat/0.33;
			stepx*=stepc;
			stepy*=stepc;
			stepz*=stepc;


		}
		else if(locrat>0.55){

			double stepc = locrat/0.55;
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
		
		coutfile << dconf.natoms << endl;
		coutfile <<"lj "<<dconf.natoms<<" U: "<<Et1<<endl;
		coutfile << setiosflags(ios::fixed) << setprecision(10);
		for(int hg = 0; hg<natoms; ++hg){
			coutfile<<"lj "<<dconf.ax[hg]<<" "<<dconf.ay[hg]<<" "<<dconf.az[hg]<<endl; 
		}
		


	}	
		

		
} //End Monte Carlo loop



		
double Uavg = Ustat.Mean();
double Usqavg = Usqstat.Mean();
double Cv = (Usqavg - Uavg*Uavg)*(Beta/Temp);
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

string itos( int Number ){
     ostringstream ss;
     ss << Number;
     return ss.str();
 }


string ftos (float input){
	ostringstream strs;
	strs << input;
	string str = strs.str();
	return str;

}

double ConfProb(double E, double Beta){

	 
	return(exp(-Beta*E));
}

__global__  void GetDoubleDeviceArrayAt_Kernel(double *d_array, int index, double * d_v){
			d_v[0]=d_array[index];

}
//get the value at index for a device vector
double GetDoubleDeviceArrayAt(double * d_array, int index){
		double h_v[1];
		double *d_v;
		h_v[0]=0.0;
		size_t sd = sizeof(double);
		cudaMalloc((void **)&d_v, sd);
		
		GetDoubleDeviceArrayAt_Kernel<<<1,1>>>(d_array, index, d_v);
		cudaDeviceSynchronize();
		cudaMemcpy(h_v, d_v, sd, cudaMemcpyDeviceToHost);
		return h_v[0];
}

__global__ void ZeroDeviceArray_Kernel(double *d_arr, int length){
	
		int ind = blockIdx.x*blockDim.x + threadIdx.x;
		
		if(ind<length){
			d_arr[ind]=0.0;
		}
}

void ZeroDeviceArray(double *d_array, int nelements){
	unsigned maxthread = MAX_THREAD;	
	//compute number of blocks needed
	unsigned nblocks = (nelements/maxthread)+1;
	ZeroDeviceArray_Kernel<<<nblocks,maxthread>>>(d_array, nelements);
	cudaDeviceSynchronize();
	return;
}

__global__ void TEloop(int natoms, int Aindex, double *box, double *ax, double *ay, double *az, double *d_Ep, int sindex){

	int bid = blockIdx.x*blockDim.x + threadIdx.x;
	
	int i = Aindex+bid+1;
	int oi = sindex+bid;

	if (i<natoms){

		double Ep=PairEnergy_Cuda(ax[i], ay[i], az[i], ax[Aindex], ay[Aindex], az[Aindex], box[0], box[1], box[2]);
		
		d_Ep[oi]=Ep;
	
		
	}	

}

double TotalEnergy_Cuda(const ConfigGPU& c){
	int maxthread = MAX_THREAD;
	double Etot = 0.0;
	int natoms = c.natoms;
	const int ti = natoms*(natoms-1)/2;
	//define and allocate device array to store pair energies
	double *d_Ep;
	size_t sti = ti*sizeof(double);
	cudaMalloc((void **)&d_Ep, sti);	
	//set all values of that array to 0.0
	ZeroDeviceArray(d_Ep, ti);
		//counter -- keep track of total number of pair interactions evaluated
		int sind = 0;
		// outer loop of pairwise interaction computation
		for(int pp = 0;pp<natoms-1;++pp){
			//number of pair interactions at this cycle of loop to compute
			int nb = natoms - pp -1;
			int nblocks = (nb/maxthread)+ (nb%maxthread == 0 ? 0 : 1);
			//parallel kernel version of inner loop		
			TEloop<<<nblocks, maxthread>>>(natoms, pp, c.d_box, c.d_ax, c.d_ay, c.d_az, d_Ep, sind);
			cudaDeviceSynchronize();
			
			sind+=nb;

		}
	//	use thrust to do parallel device prefix sum
		thrust::inclusive_scan(thrust::device, d_Ep, d_Ep+ti, d_Ep);
		//get the last value of the device energy array
		Etot=GetDoubleDeviceArrayAt(d_Ep, ti-1);
		//free up memory
		cudaFree(d_Ep);
	//return the final total potential energy
	return(Etot);
}


__global__ void DEloop(unsigned int natoms, const unsigned int Aind, double *box, double *olc, double *ax, double *ay, double *az, double *dE){

	int i = blockIdx.x*blockDim.x + threadIdx.x;
	
	
	if (i<natoms && i!=Aind){
	
		double Etot=PairEnergy_Cuda(ax[i], ay[i], az[i], ax[Aind], ay[Aind], az[Aind], box[0], box[1], box[2]);

		double Etoto=PairEnergy_Cuda(ax[i], ay[i], az[i], olc[0], olc[1], olc[2], box[0], box[1], box[2]);
		dE[i]=(double)(Etoto-Etot);
	
	}	

}

__device__ double PairEnergy_Cuda(double ax1, double ay1, double az1, double ax2, double ay2, double az2, double bx, double by, double bz){
	
	double dx, dy, dz, rs;
	double potent;
	double sigma = 1.0;
	double eps = 1.0;
	double rc = Rcut;
	double rcs = rc*rc;
	dx = ax1 - ax2;
	
	dy = ay1 - ay2;
	
	dz = az1 - az2;

	//Minimum image -- coordinates must be pre-wrapped 
	if(fabs(dx)> bx/2.0){
		dx = bx - fabs(ax1) - fabs(ax2);
	}
	
	if(fabs(dy)>by/2.0){
		dy = by - fabs(ay1) - fabs(ay2);
	}
	
	if(fabs(dz)>bz/2.0){
		dz = bz - fabs(az1) - fabs(az2);
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

// calculates the change in energy associated with one moved atom - uses gpu via cuda
double DeltaE_Cuda(const ConfigGPU& c, const unsigned& Aindex, const double& ox, const double& oy, const double& oz){
		
	int maxthread = MAX_THREAD;
	double dE=0.0;

	int natoms;
	natoms = c.natoms;
	//Allocate host copy array to store new trial coordinates
	const size_t sd3 = 3*sizeof(double);
	double * olc= (double *)malloc(sd3);

	olc[0]=ox, olc[1]=oy, olc[2]=oz;
	//Initialize and allocate the device copy	
	double * d_olc;	
	cudaMalloc((void **)&d_olc, sd3);
	// Copy input host to device
	cudaMemcpy(d_olc, olc, sd3, cudaMemcpyHostToDevice);
	//compute number of blocks needed
	int nblocks = (natoms/maxthread)+ (natoms%maxthread == 0 ? 0 : 1);
	//launch the parallel "loop" kernel
	DEloop<<<nblocks, maxthread>>>(natoms, Aindex, c.d_box, d_olc, c.d_ax, c.d_ay, c.d_az, c.d_dE);
	cudaDeviceSynchronize();
	// use thrust inclusive prefix sum on device array of dE values
	thrust::inclusive_scan(thrust::device, c.d_dE, c.d_dE+natoms, c.d_dE);
	// pull out the last element which has total
	dE=GetDoubleDeviceArrayAt(c.d_dE, natoms-1);
	// reset the d_dE array bins all to 0.0
	ZeroDeviceArray(c.d_dE, natoms);
	//free up memory 
	cudaFree(d_olc);
	free(olc);
	//return the energy difference
	return( dE ); 
	

};

