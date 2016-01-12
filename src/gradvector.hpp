//
//  gradvector.hpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/12.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//

#ifndef gradvector_hpp
#define gradvector_hpp

#include <stdio.h>

#include "field.h"
namespace nuc3d
{
    
    class gradvector
    {
        Field f_xi;
        Field f_eta;
        Field f_zeta;
        
        Field dxi;
        Field deta;
        Field dzeta;
    public:
        gradvector(int,int,int);
        void setGrad(Field &);
        ~gradvector();
    public:
        Field& getdxi(){return dxi;};
        Field& getdeta(){return deta;};
        Field& getdzeta(){return dzeta;};
        Field& getf_xi(){return f_xi;};
        Field& getf_eta(){return f_eta;};
        Field& getf_zeta(){return f_zeta;};
    };
}
#endif /* gradvector_hpp */
