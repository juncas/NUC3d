#ifndef weno5z_h
#define weno5z_h
#include "fieldOperator.h"


namespace nuc3d
{
    class weno5z : public interoplation
    {
        double ss;
        double p;
        
    public:
        weno5z();
        ~weno5z();
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
        void weno5zp(const Field & fieldIN,
                      Field & fieldOUT,
                      const int tilesize);
        
        void weno5zn(const Field & fieldIN,
                      Field & fieldOUT,
                      const int tilesize);
        
        void weno5zpBL(const Field & fieldIN,
                        const Field &boundaryL,
                        Field & fieldOUT,
                        const int tilesize);
        
        void weno5znBL(const Field & fieldIN,
                        const Field &boundaryL,
                        Field & fieldOUT,
                        const int tilesize);
        
        void weno5zpBR(const Field & fieldIN,
                        const Field &boundaryR,
                        Field & fieldOUT,
                        const int tilesize);
        
        void weno5znBR(const Field & fieldIN,
                        const Field &boundaryR,
                        Field & fieldOUT,
                        const int tilesize);


    };
}

#endif
