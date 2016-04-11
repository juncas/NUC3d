//
//  weno7js.hpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/20.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//

#ifndef weno7js_hpp
#define weno7js_hpp

#include "fieldOperator.h"


namespace nuc3d
{
    class weno7js : public interoplation
    {
        double ss;
        double p;
        
    public:
        weno7js();
        ~weno7js();
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
        void weno7jsp(const Field & fieldIN,
                      Field & fieldOUT,
                      const int tilesize);
        
        void weno7jsn(const Field & fieldIN,
                      Field & fieldOUT,
                      const int tilesize);
        
        void weno7jspBL(const Field & fieldIN,
                        const Field &boundaryL,
                        Field & fieldOUT,
                        const int tilesize);
        
        void weno7jsnBL(const Field & fieldIN,
                        const Field &boundaryL,
                        Field & fieldOUT,
                        const int tilesize);
        
        void weno7jspBR(const Field & fieldIN,
                        const Field &boundaryR,
                        Field & fieldOUT,
                        const int tilesize);
        
        void weno7jsnBR(const Field & fieldIN,
                        const Field &boundaryR,
                        Field & fieldOUT,
                        const int tilesize);
        
        
    };
}

#endif /* weno7js_hpp */
