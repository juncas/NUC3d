//
//  upwind1st.hpp
//  NUC3d
//
//  Created by Jun Peng on 15/12/12.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#ifndef upwind1st_hpp
#define upwind1st_hpp
#include "fieldOperator.h"

namespace nuc3d
{
    class upwind1st : public interoplation
    {
        
    public:
        upwind1st();
        ~upwind1st();
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

#endif /* upwind1st_hpp */
