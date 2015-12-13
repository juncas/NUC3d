#ifndef weno5z_h
#define weno5z_h
#include "fieldOperator.h"


namespace nuc3d
{
    class weno5z : public interoplation
    {
        const double coeff_weno5_alpha[3][3] = {
            {1.0/3.0,-7.0/6.0,11.0/6.0},
            {-1.0/6.0,5.0/6.0,1.0/3.0},
            {1.0/3.0,5.0/6.0,-1.0/6.0}
        };
        
        const double coeff_weno5_c[3]={0.1,0.6,0.3};
        
        const double coeff_weno5_gamma0=13.0/12.0;
        const double coeff_weno5_gamma1=1.0/4.0;
        
        double ss;
        int p;
        
    public:
        weno5z();
        ~weno5z();
        //WENO5-JS functions
        
        void interpolationInner(
                                const Field &,
                                const int,
                                const int,
                                const int,
                                const int,
                                Field &);
        void interpolationBoundaryL(
                                    const Field &,
                                    const Field &,
                                    const int,
                                    const int,
                                    const int,
                                    const int,
                                    Field &);
        void interpolationBoundaryR(
                                    const Field &,
                                    const Field &,
                                    const int,
                                    const int,
                                    const int,
                                    const int,
                                    Field &);
    private:
        double weno5zInterpolation(const double *);
    };
}

#endif
