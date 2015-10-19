#ifndef tvdrk3_h
#define tvdrk3_h
#include "fieldOperator.h"

namespace nuc3d
{
	class tvdrk3rd : public integration
	{	
        const double coeff_tvdrk3_alpha0[3][2] = { {1.0,1.0},
                                                   {3.0/4.0,1.0/4.0},
                                                   {1.0/3.0,2.0/3.0}
                                                 };
        
        std::vector<VectorField> ui;
        VectorField u;
    public:
    	tvdrk3rd();
    	~tvdrk3rd();
    	 //TVD-RK 3rd order scheme functions

        void getRHS(
                const VectorField &, // dfdx
                const VectorField &, // dgdy
                const VectorField &, // dhdz
                const double, // dt
                VectorField &);  // Right-hand-side

        void integrationAll(
                const VectorField &, // Right-hand-side: l*dt
                const VectorField &, // u(nstep)
                const VectorField &, // u_n
                int , // n th step
                VectorField &); // output: u(nstep+1)
	};
}

#endif