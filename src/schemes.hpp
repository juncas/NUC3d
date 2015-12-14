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
    double weno5jsInterpolation(const double *f,
                                 const double &ss,
                                 const double &p);
    double weno5zInterpolation(const double *f,
                                const double &ss,
                               const double &p);
    double cd6differential(const double *fl,
                           const double *fr);
    double cd2differential(const double *fl,
                           const double *fr);
    
    const double coeff_weno5_alpha[3][3] = {
        {1.0/3.0,-7.0/6.0,11.0/6.0},
        {-1.0/6.0,5.0/6.0,1.0/3.0},
        {1.0/3.0,5.0/6.0,-1.0/6.0}
    };
    
    const double coeff_weno5_c[3]={0.1,0.6,0.3};
    
    const double coeff_weno5_gamma0=13.0/12.0;
    const double coeff_weno5_gamma1=1.0/4.0;
    
    const double coeff_cd6_alpha[3] = { 3.0 / 4.0, -3.0 / 20.0,1.0 / 60.0 };

}

#endif /* schemes_hpp */
