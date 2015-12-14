#ifndef weno5js_h
#define weno5js_h
#include "fieldOperator.h"


namespace nuc3d
{
    class weno5js : public interoplation
    {        
        double ss;
        double p;
        
    public:
        weno5js();
        ~weno5js();
        //WENO5-JS functions
        
        void interpolationInner(
                                const Field &,
                                const int,
                                const int,
                                const int,
                                const int,
                                Field &,
                                const int);
        void interpolationBoundaryL(
                                    const Field &,
                                    const Field &,
                                    const int,
                                    const int,
                                    const int,
                                    const int,
                                    Field &,
                                    const int);
        void interpolationBoundaryR(
                                    const Field &,
                                    const Field &,
                                    const int,
                                    const int,
                                    const int,
                                    const int,
                                    Field &,
                                    const int);
    };
}

#endif
