//=================================
// include guard
#ifndef __ConfigGPU_H_INCLUDED__
#define __ConfigGPU_H_INCLUDED__

//=================================
// forward declared dependencies
//class Atom;
class Config;
//=================================
// included dependencies
__global__ void SetAtomCoordKernel(int index, double nx, double ny, double nz, double* d_ax, double* d_ay, double* d_az){
				d_ax[index]=nx;
				d_ay[index]=ny;
				d_az[index]=nz;
		}

__global__ void TestKernel(double* d_ax){
	int index = threadIdx.x;
	printf("atom has index %i and xcoord %f \n", index, d_ax[index]);
}


//=================================
// the actual class
class ConfigGPU {
	public:
		
		//host values
		unsigned int natoms;
		double* COM;
		double *box;
		double *ax, *ay, *az;
		size_t sz;
		size_t sd3;
		size_t si;
		double * h_dE;
		//device values
	
		//double * d_COM;
		double * d_box;
		double * d_ax, *d_ay, *d_az;
		double * d_dE;
	//Constructor
	ConfigGPU(unsigned na) {
			natoms = na;
			// Allocate host
			box = new double[3];
			COM = new double[3];
			ax = new double[na];
			ay = new double[na];
			az = new double[na];
			h_dE = new double[na];
			sz = na*sizeof(double);
			sd3 = 3*sizeof(double);
			si = sizeof(int);
			//Allocate device
			cudaMalloc((void **)&d_ax, sz);
			cudaMalloc((void **)&d_ay, sz);
			cudaMalloc((void **)&d_az, sz);
			cudaMalloc((void **)&d_box, sd3);
			cudaMalloc((void **)&d_dE, sz);
		}


		//destructor
		~ConfigGPU(){
			delete [] box;
			delete [] COM;
			delete [] ax;
			delete [] az;
			delete [] ay;
			delete [] h_dE;
			cudaFree(d_ax), cudaFree(d_ay), cudaFree(d_az);
			cudaFree(d_box);
			cudaFree(d_dE);
		}

		
				

		void FullEquateCPU(const Config& other){


			// set host copy
			double obx, oby, obz;
			//box size
			obx = (double)(other.boxx);
			oby = (double)(other.boxy);
			obz = (double)(other.boxz);
			
			box[0]=obx, box[1]=oby, box[2]=obz;
			//copy com
			double ocm;
				
			for (unsigned int i = 0; i < 3; ++i)
			{
							
				ocm = (double)(other.COM[i]);
				COM[i]=ocm;
			}
			
			
			//copy coordinates and properties
			for (unsigned int i = 0; i < natoms; ++i)
			{
				ax[i]= (double)(other.atom[i].x);
				ay[i]= (double)(other.atom[i].y);
				az[i]= (double)(other.atom[i].z);
				h_dE[i]=0.0;
			}
			
			//now make copy on device
			cudaMemcpy(d_ax, ax, sz, cudaMemcpyHostToDevice);
			cudaMemcpy(d_ay, ay, sz, cudaMemcpyHostToDevice);
			cudaMemcpy(d_az, az, sz, cudaMemcpyHostToDevice);
			cudaMemcpy(d_box, box, sd3, cudaMemcpyHostToDevice);
			cudaMemcpy(d_dE, h_dE, sz, cudaMemcpyHostToDevice);
			//cudaDeviceSynchronize();
			return;
		}
		void SetAtomCoord(int index, double nx, double ny, double nz){
				//update host copy
				ax[index]=nx, ay[index]=ny, az[index]=nz;
				//call kernel to update device copy
				SetAtomCoordKernel<<<1,1>>>(index, nx, ny, nz, d_ax, d_ay, d_az);
				cudaDeviceSynchronize();
		}

		void Test(void){

			TestKernel<<<1,1>>>(d_ax);
			cudaDeviceSynchronize();
		}

		void WrapCoord(double& xc, double& yc, double& zc){

			
			double hboxx = box[0]/2.0;
			while(xc > hboxx || xc < -hboxx) {	
				if(xc > hboxx){
					xc = xc - box[0];
				}
				else if(xc < -hboxx){
					xc = xc + box[0];
				}
			}	


	
			double hboxy = box[1]/2.0;
			while(yc > hboxy || yc < -hboxy){
				if(yc > hboxy){
					yc = yc - box[1];
				}
				else if(yc < -hboxy){
					yc = yc + box[1];
				}
			}


	
			double hboxz = box[2]/2.0;
			while(zc > hboxz || zc < -hboxz){
				if(zc > hboxz){
					zc = zc - box[2];
				}
				else if(zc < -hboxz){
					zc = zc + box[2];
				}
			}
			
			return;
	}	

};

#endif // __Config_H_INCLUDED__ 
