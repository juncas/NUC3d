//
//  schemes.hpp
//  NUC3d
//
//  Created by Jun Peng on 15/12/14.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#ifndef schemes_hpp
#define schemes_hpp

#include <stdio.h>
#include <cmath>

namespace nuc3d
{
    const double coeff_weno5_alpha[3][3] = {
        {1.0/3.0,-7.0/6.0,11.0/6.0},
        {-1.0/6.0,5.0/6.0,1.0/3.0},
        {1.0/3.0,5.0/6.0,-1.0/6.0}
    };
    
    const double coeff_weno5_c[3]={0.1,0.6,0.3};
    const double coeff_crweno5_c[3]={0.2,0.5,0.3};
    
    const double coeff_weno5_gamma0=13.0/12.0;
    const double coeff_weno5_gamma1=1.0/4.0;
        
    const double coeff_cd6_alpha[3] = { 3.0 / 4.0, -3.0 / 20.0,1.0 / 60.0 };
    
    const double coeff_cd4_alpha[2] = { 2.0 / 3.0, -1.0 / 12.0};
    
    const double coeff_cd2_alpha=0.5;

}

#endif /* schemes_hpp */
