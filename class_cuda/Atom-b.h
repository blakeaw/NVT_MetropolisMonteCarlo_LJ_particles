//=================================
// include guard
#ifndef __Atom_H_INCLUDED__
#define __Atom_H_INCLUDED__

//=================================
// forward declared dependencies


//=================================
// included dependencies


//=================================
// the actual class
class Atom {
	public:
		//int Aindex;
		double x, y, z;
		double mass;
//		Atom(Atom another) {};
		Atom(): x(0.0), y(0.0), z(0.0), mass(1.0){}
		//Copy Constructor
		Atom(const Atom& other){
			x=other.x;
			y=other.y;
			z=other.z;
			mass=other.mass;
		}
		void Equate(const Atom& other){
			x=other.x;
			y=other.y;
			z=other.z;
			mass=other.mass;
		}
};

#endif // __Atom_H_INCLUDED__ 
