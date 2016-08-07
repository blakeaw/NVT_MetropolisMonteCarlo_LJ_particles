//=================================
// include guard
#ifndef __RunningStats_H_INCLUDED__
#define __RunningStats_H_INCLUDED__

//=================================
// forward declared dependencies

//=================================
// included dependencies


//=================================
// the actual class
class RunningStats {

	public:
		RunningStats(): n(0), Mn_old(0.0), Mn_new(0.0), Sn_old(0.0), Sn_new(0.0) {}
		void Push(long double val){
			++n;
			if(n == 1){
				Mn_old = val;
				Sn_old = 0.0;
			}
			else{
				Mn_new = Mn_old + (val - Mn_old)/(double(n));
				Sn_new = Sn_old + (val - Mn_old)*(val-Mn_new);
				
				Mn_old = Mn_new;
				Sn_old = Sn_new;
				
			}
		}

	 long double Mean() const {
            return (n > 0) ? Mn_new : 0.0;
        }

        long double Variance() const {
            return ( (n > 1) ? Sn_new/(n - 1) : 0.0 );
        }

        long double StandDev() const
        {
            return sqrt( Variance() );
        }

	void Reset(void){
		n=0;
		return;
	}

    private:
        unsigned long long int n;
        long double Mn_old, Mn_new, Sn_old, Sn_new;

};


#endif // __MYCLASS_H_INCLUDED__ 
