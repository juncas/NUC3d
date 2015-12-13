#ifndef centraldifference2nd_h
#define centraldifference2nd_h
#include "fieldOperator.h"


namespace nuc3d
{
    class centraldifference2nd : public differential
    {
        const double coeff_cd6_alpha[3] = { 3.0 / 4.0, -3.0 / 20.0,1.0 / 60.0 };
        
    public:
        centraldifference2nd();
        ~centraldifference2nd();
        // 6th order central difference functions
        
        void differentialInner(
                               const Field &,
                               const int,
                               const int,
                               const int,
                               Field &);
        void differentialBoundaryL(
                                   const Field &,
                                   const Field &,
                                   const int,
                                   const int,
                                   const int,
                                   Field &);
        void differentialBoundaryR(
                                   const Field &,
                                   const Field &,
                                   const int,
                                   const int,
                                   const int,
                                   Field &);
    private:
        double cd6differential(const double *fl,
                               const double *fr);
        
    };
}



#endif