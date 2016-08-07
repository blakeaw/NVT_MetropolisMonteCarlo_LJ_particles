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
		long double x, y, z;
		long double mass;
		long double sigma;
		long double epsilon;
		string type;
//		Atom(Atom another) {};
		Atom(): x(0.0), y(0.0), z(0.0), mass(1.0), sigma(1.0), epsilon(1.0), type("LJ"){}
		//Copy Constructor
		Atom(const Atom& other){
			x=other.x;
			y=other.y;
			z=other.z;
			mass=other.mass;
			sigma=other.sigma;
			epsilon=other.epsilon;
			type=other.type;
		}
		void Equate(const Atom& other){
			x=other.x;
			y=other.y;
			z=other.z;
			mass=other.mass;
			sigma=other.sigma;
			epsilon=other.epsilon;
			type=other.type;
		}
};

#endif // __Atom_H_INCLUDED__ 
