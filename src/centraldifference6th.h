#ifndef centraldifference6th_h
#define centraldifference6th_h
#include "fieldOperator.h"


namespace nuc3d
{
    class centraldifference6th : public differential
    {
        
    public:
        centraldifference6th();
        ~centraldifference6th();
        // 6th order central difference functions
        
        void differentialInner(const Field &,
                               Field &,
                               const int);
        
        void differentialBoundaryL(const Field &,
                                   const Field &,
                                   Field &,
                                   const int);
        void differentialBoundaryR(const Field &,
                                   const Field &,
                                   Field &,
                                   const int);        
    };
}



#endif