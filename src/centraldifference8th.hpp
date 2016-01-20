//
//  centraldifferebce8th.hpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/20.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//

#ifndef centraldifference8th_hpp
#define centraldifference8th_hpp

#include "fieldOperator.h"


namespace nuc3d
{
    class centraldifference8th : public differential
    {
        
    public:
        centraldifference8th();
        ~centraldifference8th();
        // 8th order central difference functions
        
        void differentialInner(const Field &,
                               Field &,
                               const int);
        
        void differentialBoundaryL(const Field &,
                                   const Field &,
                                   Field &,
                                   const int);
        void differentialBoundaryR(const Field &,
                                   const Field &,
                                   Field &,
                                   const int);
    };
}


#endif /* centraldifference8th_hpp */
