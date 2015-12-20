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
        void upwind1stp(const Field & fieldIN,
                      Field & fieldOUT,
                      const int tilesize);
        
        void upwind1stn(const Field & fieldIN,
                      Field & fieldOUT,
                      const int tilesize);
        
        void upwind1stpBL(const Field & fieldIN,
                        const Field &boundaryL,
                        Field & fieldOUT,
                        const int tilesize);
        
        void upwind1stnBL(const Field & fieldIN,
                        const Field &boundaryL,
                        Field & fieldOUT,
                        const int tilesize);
        
        void upwind1stpBR(const Field & fieldIN,
                        const Field &boundaryR,
                        Field & fieldOUT,
                        const int tilesize);
        
        void upwind1stnBR(const Field & fieldIN,
                        const Field &boundaryR,
                        Field & fieldOUT,
                        const int tilesize);


    };
}

#endif /* upwind1st_hpp */
