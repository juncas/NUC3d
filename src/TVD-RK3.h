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
        
        VectorField u0;
        std::vector<VectorField> ui;
    public:
        tvdrk3rd();
        ~tvdrk3rd();
        //TVD-RK 3rd order scheme functions
        
        void initial(const VectorField &);
        
        void integrationAll(const VectorField &, // Right-hand-side: l*dt
                            VectorField &, // u_n
                            double,
                            int); // step n
    private:
        
        void rk1st(const VectorField &,
                   VectorField &,
                   double );
        void rk2nd(const VectorField &,
                   VectorField &,
                   double );
        void rk3rd(const VectorField &,
                   VectorField &,
                   double );
    };
}

#endif