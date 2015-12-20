//
//  schemes.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/12/14.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#include "schemes.hpp"

void nuc3d::weno5jsInterpolation(double &rf,
                                   const double *f,const double &ss,const double &p)
{
    double omega0,omega1,omega2;
    double alpha0,alpha1,alpha2,alphaSum;
    double is0,is1,is2;
    double q30,q31,q32;
    
    is0= coeff_weno5_gamma0*pow((    f[0]-2.0*f[1]+    f[2]),2)
    +coeff_weno5_gamma1*pow((    f[0]-4.0*f[1]+3.0*f[2]),2);
    
    is1= coeff_weno5_gamma0*pow((    f[1]-2.0*f[2]+    f[3]),2)
    +coeff_weno5_gamma1*pow((    f[1]-             f[3]),2);
    
    is2= coeff_weno5_gamma0*pow((    f[2]-2.0*f[3]+    f[4]),2)
    +coeff_weno5_gamma1*pow((3.0*f[2]-4.0*f[3]+    f[4]),2);
    
    
    q30= coeff_weno5_alpha[0][0]*f[0]
    +coeff_weno5_alpha[0][1]*f[1]
    +coeff_weno5_alpha[0][2]*f[2];
    
    q31= coeff_weno5_alpha[1][0]*f[1]
    +coeff_weno5_alpha[1][1]*f[2]
    +coeff_weno5_alpha[1][2]*f[3];
    
    q32= coeff_weno5_alpha[2][0]*f[2]
    +coeff_weno5_alpha[2][1]*f[3]
    +coeff_weno5_alpha[2][2]*f[4];
    
    alpha0=coeff_weno5_c[0]/pow((ss+is0),p);
    alpha1=coeff_weno5_c[1]/pow((ss+is1),p);
    alpha2=coeff_weno5_c[2]/pow((ss+is2),p);
    
    
    alphaSum=alpha0+alpha1+alpha2;
    
    omega0=alpha0/alphaSum;
    omega1=alpha1/alphaSum;
    omega2=alpha2/alphaSum;
    
    rf=omega0*q30+omega1*q31+omega2*q32;
}

double nuc3d::weno5zInterpolation(const double *f,const double &ss,const double &p)

{
    double omega0,omega1,omega2;
    double alpha0,alpha1,alpha2,alphaSum;
    double is0,is1,is2;
    double q30,q31,q32;
    double tau;
    
    is0= coeff_weno5_gamma0*pow((    f[0]-2.0*f[1]+    f[2]),2)
    +coeff_weno5_gamma1*pow((    f[0]-4.0*f[1]+3.0*f[2]),2);
    
    is1= coeff_weno5_gamma0*pow((    f[1]-2.0*f[2]+    f[3]),2)
    +coeff_weno5_gamma1*pow((    f[1]-             f[3]),2);
    
    is2= coeff_weno5_gamma0*pow((    f[2]-2.0*f[3]+    f[4]),2)
    +coeff_weno5_gamma1*pow((3.0*f[2]-4.0*f[3]+    f[4]),2);
    
    
    q30= coeff_weno5_alpha[0][0]*f[0]
    +coeff_weno5_alpha[0][1]*f[1]
    +coeff_weno5_alpha[0][2]*f[2];
    
    q31= coeff_weno5_alpha[1][0]*f[1]
    +coeff_weno5_alpha[1][1]*f[2]
    +coeff_weno5_alpha[1][2]*f[3];
    
    q32= coeff_weno5_alpha[2][0]*f[2]
    +coeff_weno5_alpha[2][1]*f[3]
    +coeff_weno5_alpha[2][2]*f[4];
    
    
    tau=std::abs(is2-is0);
    
    alpha0=coeff_weno5_c[0]*(1.0+std::pow(tau/(is0+ss),p));
    alpha1=coeff_weno5_c[1]*(1.0+std::pow(tau/(is1+ss),p));
    alpha2=coeff_weno5_c[2]*(1.0+std::pow(tau/(is2+ss),p));
    
    alphaSum=alpha0+alpha1+alpha2;
    
    omega0=alpha0/alphaSum;
    omega1=alpha1/alphaSum;
    omega2=alpha2/alphaSum;
    
    return omega0*q30+omega1*q31+omega2*q32;
}


double nuc3d::cd6differential(const double *fl,
                              const double *fr)
{
    double df;
    
    df = coeff_cd6_alpha[2]*(fr[2]-fl[2])
    +coeff_cd6_alpha[1]*(fr[1]-fl[1])
    +coeff_cd6_alpha[0]*(fr[0]-fl[0]);
    
    return df;
}


double nuc3d::cd2differential(const double &fl,
                              const double &fr)
{
    double df;
    
    df = 0.5*(fr-fl);
    
    return df;
}
