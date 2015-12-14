#ifndef weno5z_h
#define weno5z_h
#include "fieldOperator.h"


namespace nuc3d
{
    class weno5z : public interoplation
    {        
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
                                Field &,
                                const int tilesize);
        void interpolationBoundaryL(
                                    const Field &,
                                    const Field &,
                                    const int,
                                    const int,
                                    const int,
                                    const int,
                                    Field &,
                                    const int tilesize);
        void interpolationBoundaryR(
                                    const Field &,
                                    const Field &,
                                    const int,
                                    const int,
                                    const int,
                                    const int,
                                    Field &,
                                    const int tilesize);

    };
}

#endif
