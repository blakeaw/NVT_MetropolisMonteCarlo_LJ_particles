//=================================
// include guard
#ifndef __Config_H_INCLUDED__
#define __Config_H_INCLUDED__

//=================================
// forward declared dependencies
class Atom;
//=================================
// included dependencies


//=================================
// the actual class
class Config {
	public:
		
		unsigned int natoms;
		double* COM;
		double boxx, boxy, boxz;
		double boxold;
		Atom * atom;
		
		//Config(): natoms(a_natoms_a), boxx(1.0), boxy(1.0), boxz(1.0) {
 		//	atom = new Atom[a_natoms_a];
		//	COM = new double[3];
		//}

		Config(unsigned na): boxx(0.0), boxy(0.0), boxz(0.0){
			natoms = na;
			//boxx = 0.0, boxy = 0.0, boxz = 0.0;
			atom = new Atom[na];
			COM = new double[3];
	
		}

		Config(void): boxx(0.0), boxy(0.0), boxz(0.0){
			//boxx = 0.0, boxy = 0.0, boxz = 0.0;
			
			COM = new double[3];
	
		}

		~Config(){
			delete [] atom;
			delete [] COM;
		}

		//Copy Constructor

		Config(const Config& other){
			natoms = other.natoms;
			atom = new Atom[natoms];
			COM = new double[3];

			//box size
			boxx = other.boxx;
			boxy = other.boxy;
			boxz = other.boxz;
			
			//copy com
			for (unsigned int i = 0; i < 3; ++i)
			{
				COM[i] = other.COM[i];
			}
			//copy coordinates and properties
			for (unsigned int i = 0; i < natoms; ++i)
			{
				//call atom equate
				atom[i].Equate(other.atom[i]);
			}
			
			return;
		}
				
		void Equate(const Config& other){
			
			for(unsigned int i = 0; i<natoms;++i){
				atom[i].x = other.atom[i].x;
				atom[i].y = other.atom[i].y;
				atom[i].z = other.atom[i].z;
			}
		}
		
		void FullEquate(const Config& other){
			
			//box size
			boxx = other.boxx;
			boxy = other.boxy;
			boxz = other.boxz;
			
			//copy com
			for (unsigned int i = 0; i < 3; ++i)
			{
				COM[i] = other.COM[i];
			}
			//copy coordinates and properties
			for (unsigned int i = 0; i < natoms; ++i)
			{
				//call atom equate
				atom[i].Equate(other.atom[i]);
			}
			
			return;
		}

		void CalcCom(void){
			
			double M = 0.0;
			double Rx, Ry, Rz;
			double sumrx=0.0, sumry=0.0, sumrz=0.0;
			
			
			for(unsigned int i = 0;i<natoms;++i){
				
				double mass = atom[i].mass;
				sumrx += atom[i].x*mass;
				sumry += atom[i].y*mass;
				sumrz += atom[i].z*mass;
				M += mass;
				
			}
			Rx = sumrx/M;
			Ry = sumry/M;
			Rz = sumrz/M;

			COM[0]=Rx;
			COM[1]=Ry;
			COM[2]=Rz;
			return;

		}
	
		double RadGyr(void){
			CalcCom();
			double sumgyrx = 0.0, sumgyry = 0.0, sumgyrz = 0.0;
			double M = 0.0;
			double Rsq, R;
			double sumgyr=0.0;
//			cout << " M " << M << endl;
			for(unsigned long y = 0;y<natoms; ++y){
				double rx, ry, rz;
				rx = atom[y].x;
				ry = atom[y].y;
				rz = atom[y].z;
				double mass = atom[y].mass;
				sumgyrx = (rx - COM[0])*(rx - COM[0]);
				sumgyry = (ry - COM[1])*(ry - COM[1]);
				sumgyrz = (rz - COM[2])*(rz - COM[2]);
				sumgyr += (sumgyrx + sumgyry + sumgyrz)*mass;
				M+= mass;		
			}

			Rsq = sumgyr;
		
			R = sqrt(Rsq/M);
	
			return(R);	
		}
		
		void ScaleCoordinates(double& Vscale){

			
			double Cscale = pow(Vscale, 1.0/3.0);
			for (unsigned i = 0; i < natoms; ++i)
			{
				atom[i].x *= Cscale;
				atom[i].y *= Cscale;
				atom[i].z *= Cscale;
			}


			return;

		}
		void ScaleCoordinates(double xscale, double yscale, double zscale){

			for (unsigned i = 0; i < natoms; ++i)
			{
				atom[i].x *= xscale;
				atom[i].y *= yscale;
				atom[i].z *= zscale;
			}


			return;

		}
		double ComDist(unsigned Aindex){
			double dx, dy, dz;
			dx = atom[Aindex].x - COM[0];
			dy = atom[Aindex].y - COM[0];
			dz = atom[Aindex].z - COM[0];
			double r, rsq;
			rsq = (dx*dx) + (dy*dy) + (dz*dz);
			r = sqrt(rsq);

			return(r);
}

		bool CheckComDistAll(const double& dist){
			CalcCom();
			bool flag = 1;
			for (unsigned int i = 0; i < natoms; ++i)
			{
				double dx, dy, dz;
				dx = atom[i].x - COM[0];
				dy = atom[i].y - COM[0];
				dz = atom[i].z - COM[0];
				double r, rsq;
				rsq = (dx*dx) + (dy*dy) + (dz*dz);
				r = sqrt(rsq);
				if (r>dist){
					flag=0;
					break;
				}
			}
			

			return flag;
}

	void SetBox(double& bx, double& by, double& bz){

		boxold=boxx;
		boxx=bx;
		boxy=by;
		boxz=bz;	
		return;
	}
	
	void InitializeAtoms(unsigned nat){
		natoms = nat;
		atom = new Atom[nat];
		
	}

	void WrapCoordinates(void){

		for (unsigned int i = 0; i < natoms; ++i)
		{
			double xc = atom[i].x;
			double yc = atom[i].y;
			double zc = atom[i].z;

			double hboxx = boxx/2.0;
			while(xc > hboxx || xc < -hboxx) {	
				if(xc > hboxx){
					xc = xc - boxx;
				}
				else if(xc < -hboxx){
					xc = xc + boxx;
				}
			}	


	
			double hboxy = boxy/2.0;
			while(yc > hboxy || yc < -hboxy){
				if(yc > hboxy){
					yc = yc - boxy;
				}
				else if(yc < -hboxy){
					yc = yc + boxy;
				}
			}


	
			double hboxz = boxz/2.0;
			while(zc > hboxz || zc < -hboxz){
				if(zc > hboxz){
					zc = zc - boxz;
				}
				else if(zc < -hboxz){
					zc = zc + boxz;
				}
			}

			atom[i].x = xc;
			atom[i].y = yc;
			atom[i].z = zc;
		}
		
		return;

}	

	void WrapAtom(const unsigned& aindex){

			double xc = atom[aindex].x;
			double yc = atom[aindex].y;
			double zc = atom[aindex].z;

			double hboxx = boxx/2.0;
			while(xc > hboxx || xc < -hboxx) {	
				if(xc > hboxx){
					xc = xc - boxx;
				}
				else if(xc < -hboxx){
					xc = xc + boxx;
				}
			}	


	
			double hboxy = boxy/2.0;
			while(yc > hboxy || yc < -hboxy){
				if(yc > hboxy){
					yc = yc - boxy;
				}
				else if(yc < -hboxy){
					yc = yc + boxy;
				}
			}


	
			double hboxz = boxz/2.0;
			while(zc > hboxz || zc < -hboxz){
				if(zc > hboxz){
					zc = zc - boxz;
				}
				else if(zc < -hboxz){
					zc = zc + boxz;
				}
			}

			atom[aindex].x = xc;
			atom[aindex].y = yc;
			atom[aindex].z = zc;
			
			return;
	}	
	
	void WrapCoord(double& xc, double& yc, double& zc){

			
			double hboxx = boxx/2.0;
			while(xc > hboxx || xc < -hboxx) {	
				if(xc > hboxx){
					xc = xc - boxx;
				}
				else if(xc < -hboxx){
					xc = xc + boxx;
				}
			}	


	
			double hboxy = boxy/2.0;
			while(yc > hboxy || yc < -hboxy){
				if(yc > hboxy){
					yc = yc - boxy;
				}
				else if(yc < -hboxy){
					yc = yc + boxy;
				}
			}


	
			double hboxz = boxz/2.0;
			while(zc > hboxz || zc < -hboxz){
				if(zc > hboxz){
					zc = zc - boxz;
				}
				else if(zc < -hboxz){
					zc = zc + boxz;
				}
			}
			
			return;
	}	
	
	void ReCenterZero(void){
		for (unsigned int i = 0; i < natoms; ++i)
		{
			atom[i].x += -COM[0];
			atom[i].y += -COM[1];
			atom[i].z +=-COM[2];
		}
	}
};

#endif // __Config_H_INCLUDED__ 
