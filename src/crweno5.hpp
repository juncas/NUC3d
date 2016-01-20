//
//  crweno5.hpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/20.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//

#ifndef crweno5_hpp
#define crweno5_hpp
#include "fieldOperator.h"

namespace nuc3d
{
    class crweno5 : public interoplation
    {
        double ss;
        double p;
        
    public:
        crweno5();
        ~crweno5();
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
        void crweno5p(const Field & fieldIN,
                      Field & fieldOUT,
                      const int tilesize);
        
        void crweno5n(const Field & fieldIN,
                      Field & fieldOUT,
                      const int tilesize);
        
        void crweno5pBL(const Field & fieldIN,
                        const Field &boundaryL,
                        Field & fieldOUT,
                        const int tilesize);
        
        void crweno5nBL(const Field & fieldIN,
                        const Field &boundaryL,
                        Field & fieldOUT,
                        const int tilesize);
        
        void crweno5pBR(const Field & fieldIN,
                        const Field &boundaryR,
                        Field & fieldOUT,
                        const int tilesize);
        
        void crweno5nBR(const Field & fieldIN,
                        const Field &boundaryR,
                        Field & fieldOUT,
                        const int tilesize);
        
        
    };
}


#endif /* crweno5_hpp */
