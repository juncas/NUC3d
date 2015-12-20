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
        
        void interpolationInner(const Field &,
                                const int,
                                Field &,
                                const int);
        
        void interpolationBoundaryL(const Field &,
                                    const Field &,
                                    const int,
                                    Field &,
                                    const int);
        
        void interpolationBoundaryR(const Field &,
                                    const Field &,
                                    const int,
                                    Field &,
                                    const int);
    private:
        void weno5jsp(const Field & fieldIN,
                      Field & fieldOUT,
                      const int tilesize);
        
        void weno5jsn(const Field & fieldIN,
                      Field & fieldOUT,
                      const int tilesize);
        
        void weno5jspBL(const Field & fieldIN,
                        const Field &boundaryL,
                        Field & fieldOUT,
                        const int tilesize);
        
        void weno5jsnBL(const Field & fieldIN,
                        const Field &boundaryL,
                        Field & fieldOUT,
                        const int tilesize);
        
        void weno5jspBR(const Field & fieldIN,
                        const Field &boundaryR,
                        Field & fieldOUT,
                        const int tilesize);
        
        void weno5jsnBR(const Field & fieldIN,
                        const Field &boundaryR,
                        Field & fieldOUT,
                        const int tilesize);
        
        
    };
}

#endif
