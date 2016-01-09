//
//  boundaryScheme.hpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/9.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//

#ifndef boundaryScheme_hpp
#define boundaryScheme_hpp

#include <stdio.h>

#include "fieldOperator.h"


namespace nuc3d
{
    class boundaryScheme: public differential
    {
        typedef void (boundaryScheme::*pBCscheme)(double *,double *,int &);
    public:
        pBCscheme myBCschemeL[4]={&boundaryScheme::firstOrderL,
            &boundaryScheme::secondOrder,
            &boundaryScheme::fourthOrder,
            &boundaryScheme::sixthOrder
        };
        
        pBCscheme myBCschemeR[4]={&boundaryScheme::firstOrderR,
            &boundaryScheme::secondOrder,
            &boundaryScheme::fourthOrder,
            &boundaryScheme::sixthOrder
        };
        boundaryScheme();
        ~boundaryScheme();
        // 6th order central difference functions
        
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
    private:
        void firstOrderL(double *,double *,int &);
        void firstOrderR(double *,double *,int &);
        void secondOrder(double *,double *,int &);
        void fourthOrder(double *,double *,int &);
        void sixthOrder(double *,double *,int &);
    };
}




#endif /* boundaryScheme_hpp */
