//=================================
// include guard
#ifndef __MTRNG_H_INCLUDED__
#define __MTRNG_H_INCLUDED__

//=================================
// forward declared dependencies

//=================================
// included dependencies
// Values for RNG
#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL /* Least significant 31 bits */

//=================================
// the actual class
class MTRandomNum {
      public:
		
	unsigned long long int mt[NN];
    int mti;


	// initializes with a seed
    void initialize(unsigned long long int seed){
        	mt[0] = seed;
   		 for (mti=1; mti<NN; mti++){ 
      			  mt[mti] =  (unsigned long long int)(6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
		}
	
		return;
        }
	// called to generate random number
	// double on [0,1] 
	long double generate(void){
		 int i;
   		 unsigned long long x;
   		 static unsigned long long mag01[2]={0ULL, MATRIX_A};

    		if (mti >= NN) { /* generate NN words at one time */	

        	/* if initialize has not been called, */
        	/* a default initial seed is used     */
        		if (mti == NN+1){
        		    initialize(5489ULL);
			} 

        		for (i=0;i<NN-MM;i++) {
        		    x = (mt[i]&UM)|(mt[i+1]&LM);
        		    mt[i] = mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        		}
        		for (;i<NN-1;i++) {
        		    x = (mt[i]&UM)|(mt[i+1]&LM);
        		    mt[i] = mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        		}
        		x = (mt[NN-1]&UM)|(mt[0]&LM);
        		mt[NN-1] = mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
	
        		mti = 0;
    		}
  
   		 x = mt[mti++];

    		x ^= (x >> 29) & 0x5555555555555555ULL;
    		x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    		x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    		x ^= (x >> 43);
		// [0,1] real
//		return (x >> 11) * (1.0/9007199244640991.0);
		// [0,1) real
//		return (x >> 11) * (1.0/9007199254740992.0);
		// (0,1) real 
		return ((x >> 12) + 0.5) * (1.0/4503599627370496.0);
	
	}
};


#endif // __MTRNG_H_INCLUDED__ 
