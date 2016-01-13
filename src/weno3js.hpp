//
//  weno3js.hpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/13.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//

#ifndef weno3js_hpp
#define weno3js_hpp

#include "fieldOperator.h"


namespace nuc3d
{
    class weno3js : public interoplation
    {
        double ss;
        int p;
        
    public:
        weno3js();
        ~weno3js();
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
        void weno3jsp(const Field & fieldIN,
                      Field & fieldOUT,
                      const int tilesize);
        
        void weno3jsn(const Field & fieldIN,
                      Field & fieldOUT,
                      const int tilesize);
        
        void weno3jspBL(const Field & fieldIN,
                        const Field &boundaryL,
                        Field & fieldOUT,
                        const int tilesize);
        
        void weno3jsnBL(const Field & fieldIN,
                        const Field &boundaryL,
                        Field & fieldOUT,
                        const int tilesize);
        
        void weno3jspBR(const Field & fieldIN,
                        const Field &boundaryR,
                        Field & fieldOUT,
                        const int tilesize);
        
        void weno3jsnBR(const Field & fieldIN,
                        const Field &boundaryR,
                        Field & fieldOUT,
                        const int tilesize);
        
        
    };
}
#endif /* weno3js_hpp */
