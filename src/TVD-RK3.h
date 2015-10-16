#ifndef tvdrk3_h
#define tvdrk3_h
#include "fieldOperator.h"

namespace nuc3d
{
	class tvdrk3 : public integration
	{	
        const double coeff_tvdrk3_alpha0[3][2] = { {1.0,1.0},
                                                   {3.0/4.0,1.0/4.0},
                                                   {1.0/3.0,2.0/3.0}
                                                 };

    public:
    	tvdrk3();
    	~tvdrk3();
    	 //TVD-RK 3rd order scheme functions

        void getRHS(
                const Field &, // dfdx
                const Field &, // dgdy
                const Field &, // dhdz
                const double, // dt
                Field &);  // Right-hand-side

        void integrationAll(
                const Field &, // Right-hand-side: l*dt
                const Field &, // u(nstep)
                const Field &, // u_n
                int , // n th step
                Field &); // output: u(nstep+1)
	};
}

#endif