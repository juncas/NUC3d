//
//  hccs.hpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/24.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//

#ifndef hccs_hpp
#define hccs_hpp
#include "fieldOperator.h"

namespace nuc3d
{
    class hccs : public interoplation
    {
        double ss;
        double p;
        
    public:
        hccs();
        ~hccs();
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
        void hccsp(const Field & fieldIN,
                      Field & fieldOUT,
                      const int tilesize);
        
        void hccsn(const Field & fieldIN,
                      Field & fieldOUT,
                      const int tilesize);
        
        void hccspBL(const Field & fieldIN,
                        const Field &boundaryL,
                        Field & fieldOUT,
                        const int tilesize);
        
        void hccsnBL(const Field & fieldIN,
                        const Field &boundaryL,
                        Field & fieldOUT,
                        const int tilesize);
        
        void hccspBR(const Field & fieldIN,
                        const Field &boundaryR,
                        Field & fieldOUT,
                        const int tilesize);
        
        void hccsnBR(const Field & fieldIN,
                        const Field &boundaryR,
                        Field & fieldOUT,
                        const int tilesize);
        
        
    };
}

#endif /* hccs_hpp */
