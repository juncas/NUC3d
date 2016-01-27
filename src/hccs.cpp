//
//  hccs.cpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/24.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//

#include "hccs.hpp"

#include "schemes.hpp"

nuc3d::hccs::hccs():
ss(1.0e-40),
p(2.0)
{
    
}

nuc3d::hccs::~hccs()
{
    
}

void nuc3d::hccs::interpolationInner(const Field & fieldIN,
                                     const int uw,//uw=(1,-1)
                                     Field & fieldOUT,
                                     const int tilesize)
{
    switch(uw)
    {
        case 1:
            hccsp(fieldIN, fieldOUT, tilesize);
            break;
        case -1:
            hccsn(fieldIN, fieldOUT, tilesize);
            break;
        default:
            std::cout<<"weno error: no such direction"<<std::endl;
            exit(-1);
    }
    
}

void nuc3d::hccs::hccsp(const Field & fieldIN,
                        Field & fieldOUT,
                        const int tilesize)
{
    
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    
    int nx=fieldIN.getSizeX();
    int ny=fieldIN.getSizeY();
    int nz=fieldIN.getSizeZ();
    
    int nx0=fieldOUT.getSizeX();
    int ny0=fieldOUT.getSizeY();
    int nz0=fieldOUT.getSizeZ();
    //    vector is very very slow as this function is called many times
    //    std::vector<double> a(nx-2*tilesize+1,0.0);
    //    std::vector<double> b(nx-2*tilesize+1,1.0);
    //    std::vector<double> c(nx-2*tilesize+1,0.0);
    //    std::vector<double> d(nx-2*tilesize+1,0.0);
    double a[256];
    double b[256];
    double c[256];
    double d[256];
    
    
    int ibeg=(tilesize-1);
    int iend=nx-tilesize;
    int jbeg=0;
    int jend=ny;
    int kbeg=0;
    int kend=nz;
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<(ibeg+1);i++)
            {
                double omega0,omega1,omega2,omega3;
                double alpha0,alpha1,alpha2,alpha3,alphaSum;
                double is0,is1,is2,is3;
                double q40,q41,q42,q43;
                double tau7;
                
                int idx_f=nx*ny*k+nx*j+i;
                
                is0=pIn[idx_f-3]*(547.0*pIn[idx_f-3]-3882.0*pIn[idx_f-2]+4642.0*pIn[idx_f-1]-1854.0*pIn[idx_f])
                +pIn[idx_f-2]*(7043.0*pIn[idx_f-2]-17246.0*pIn[idx_f-1]+7042.0*pIn[idx_f])
                +pIn[idx_f-1]*(11003.0*pIn[idx_f-1]-9402.0*pIn[idx_f])
                +2107.0*pIn[idx_f]*pIn[idx_f];
                
                is1=pIn[idx_f-2]*(267.0*pIn[idx_f-2]-1642.0*pIn[idx_f-1]+1602.0*pIn[idx_f]-494.0*pIn[idx_f+1])
                +pIn[idx_f-1]*(2843.0*pIn[idx_f-1]-5966.0*pIn[idx_f]+1922.0*pIn[idx_f+1])
                +pIn[idx_f]*(3443.0*pIn[idx_f]-2522.0*pIn[idx_f+1])
                +547.0*pIn[idx_f+1]*pIn[idx_f+1];
                
                
                is2=pIn[idx_f-1]*(547.0*pIn[idx_f-1]-2522.0*pIn[idx_f]+1922.0*pIn[idx_f+1]-494.0*pIn[idx_f+2])
                +pIn[idx_f]*(3443.0*pIn[idx_f]-5966.0*pIn[idx_f+1]+1602.0*pIn[idx_f+2])
                +pIn[idx_f+1]*(2843.0*pIn[idx_f+1]-1642.0*pIn[idx_f+2])
                +267.0*pIn[idx_f+2]*pIn[idx_f+2];
                
                is3=pIn[idx_f]*(2107.0*pIn[idx_f]-9402.0*pIn[idx_f+1]+7042.0*pIn[idx_f+2]-1854.0*pIn[idx_f+3])
                +pIn[idx_f+1]*(11003.0*pIn[idx_f+1]-17246.0*pIn[idx_f+2]+4642.0*pIn[idx_f+3])
                +pIn[idx_f+2]*(7043.0*pIn[idx_f+2]-3882.0*pIn[idx_f+3])
                +547.0*pIn[idx_f+3]*pIn[idx_f+3];
                
                q40=coeff_weno7_alpha[0][0]*pIn[idx_f-3]+coeff_weno7_alpha[0][1]*pIn[idx_f-2]
                +coeff_weno7_alpha[0][2]*pIn[idx_f-1]+coeff_weno7_alpha[0][3]*pIn[idx_f];
                
                q41=coeff_weno7_alpha[1][0]*pIn[idx_f-2]+coeff_weno7_alpha[1][1]*pIn[idx_f-1]
                +coeff_weno7_alpha[1][2]*pIn[idx_f]+coeff_weno7_alpha[1][3]*pIn[idx_f+1];
                
                q42=coeff_weno7_alpha[2][0]*pIn[idx_f-1]+coeff_weno7_alpha[2][1]*pIn[idx_f]
                +coeff_weno7_alpha[2][2]*pIn[idx_f+1]+coeff_weno7_alpha[2][3]*pIn[idx_f+2];
                
                q43=coeff_weno7_alpha[3][0]*pIn[idx_f]+coeff_weno7_alpha[3][1]*pIn[idx_f+1]
                +coeff_weno7_alpha[3][2]*pIn[idx_f+2]+coeff_weno7_alpha[3][3]*pIn[idx_f+3];
                
                tau7=std::abs(is0-is3);
                alpha0=coeff_weno7_c[0]*(1.0+std::pow(tau7/(ss+is0),2.0));
                alpha1=coeff_weno7_c[1]*(1.0+std::pow(tau7/(ss+is1),2.0));
                alpha2=coeff_weno7_c[2]*(1.0+std::pow(tau7/(ss+is2),2.0));
                alpha3=coeff_weno7_c[3]*(1.0+std::pow(tau7/(ss+is3),2.0));
                
                
                alphaSum=alpha0+alpha1+alpha2+alpha3;
                
                omega0=alpha0/alphaSum;
                omega1=alpha1/alphaSum;
                omega2=alpha2/alphaSum;
                omega3=alpha3/alphaSum;
                
                a[i-ibeg]=0.0;
                b[i-ibeg]=1.0;
                c[i-ibeg]=0.0;
                d[i-ibeg]=omega0*q40+omega1*q41+omega2*q42+omega3*q43;
            }
            
            for(int i=ibeg+1;i<(iend-1);i++)
            {
                double omega0,omega1,omega2;
                double alpha0,alpha1,alpha2,alphaSum;
                double is0,is1,is2;
                double q7;
                
                int idx_f=nx*ny*k+nx*j+i;
                
                is0= coeff_weno5_gamma0*std::pow((pIn[idx_f-2]-2.0*pIn[idx_f-1]+pIn[idx_f]),2)
                +coeff_weno5_gamma1*std::pow((pIn[idx_f-2]-4.0*pIn[idx_f-1]+3.0*pIn[idx_f]),2);
                
                is1= coeff_weno5_gamma0*std::pow((pIn[idx_f-1]-2.0*pIn[idx_f]+pIn[idx_f+1]),2)
                +coeff_weno5_gamma1*std::pow((pIn[idx_f-1]-pIn[idx_f+1]),2);
                
                is2= coeff_weno5_gamma0*std::pow((pIn[idx_f]-2.0*pIn[idx_f+1]+pIn[idx_f+2]),2)
                +coeff_weno5_gamma1*std::pow((3.0*pIn[idx_f]-4.0*pIn[idx_f+1]+pIn[idx_f+2]),2);
                
                q7=(-pIn[idx_f-2]+19.0*pIn[idx_f-1]+239.0*pIn[idx_f]+159.0*pIn[idx_f+1]+4.0*pIn[idx_f+2])/420.0;
                
                double tau=std::abs(is2-is0);
                
                alpha0=coeff_crweno5_c[0]*(1.0+std::pow(tau/(is0+ss),p));
                alpha1=coeff_crweno5_c[1]*(1.0+std::pow(tau/(is1+ss),p));
                alpha2=coeff_crweno5_c[2]*(1.0+std::pow(tau/(is2+ss),p));
                
                alphaSum=alpha0+alpha1+alpha2;
                
                omega0=alpha0/alphaSum;
                omega1=alpha1/alphaSum;
                omega2=alpha2/alphaSum;
                
                double theta=1.0/(1.0+std::pow(alphaSum-1.0,2.0));
                
                a[i-ibeg]=(2.0/3.0*omega0+1.0/3.0*omega1)*(1.0-theta)
                +2.0/7.0*theta;
                b[i-ibeg]=(1.0/3.0*omega0+2.0/3.0*(omega1+omega2))*(1.0-theta)
                +4.0/7.0*theta;
                c[i-ibeg]=(1.0/3.0*omega2)*(1.0-theta)
                +1.0/7.0*theta;
                d[i-ibeg]=(1.0/6.0*omega0*pIn[idx_f-1]
                           +(5.0*(omega0+omega1)+omega2)*pIn[idx_f]/6.0
                           +(omega1+5.0*omega2)*pIn[idx_f+1]/6.0)*(1.0-theta)+q7*theta;
            }
            
            for(int i=iend-1;i<iend;i++)
            {
                double omega0,omega1,omega2,omega3;
                double alpha0,alpha1,alpha2,alpha3,alphaSum;
                double is0,is1,is2,is3;
                double q40,q41,q42,q43;
                double tau7;
                
                int idx_f=nx*ny*k+nx*j+i;
                
                is0=pIn[idx_f-3]*(547.0*pIn[idx_f-3]-3882.0*pIn[idx_f-2]+4642.0*pIn[idx_f-1]-1854.0*pIn[idx_f])
                +pIn[idx_f-2]*(7043.0*pIn[idx_f-2]-17246.0*pIn[idx_f-1]+7042.0*pIn[idx_f])
                +pIn[idx_f-1]*(11003.0*pIn[idx_f-1]-9402.0*pIn[idx_f])
                +2107.0*pIn[idx_f]*pIn[idx_f];
                
                is1=pIn[idx_f-2]*(267.0*pIn[idx_f-2]-1642.0*pIn[idx_f-1]+1602.0*pIn[idx_f]-494.0*pIn[idx_f+1])
                +pIn[idx_f-1]*(2843.0*pIn[idx_f-1]-5966.0*pIn[idx_f]+1922.0*pIn[idx_f+1])
                +pIn[idx_f]*(3443.0*pIn[idx_f]-2522.0*pIn[idx_f+1])
                +547.0*pIn[idx_f+1]*pIn[idx_f+1];
                
                
                is2=pIn[idx_f-1]*(547.0*pIn[idx_f-1]-2522.0*pIn[idx_f]+1922.0*pIn[idx_f+1]-494.0*pIn[idx_f+2])
                +pIn[idx_f]*(3443.0*pIn[idx_f]-5966.0*pIn[idx_f+1]+1602.0*pIn[idx_f+2])
                +pIn[idx_f+1]*(2843.0*pIn[idx_f+1]-1642.0*pIn[idx_f+2])
                +267.0*pIn[idx_f+2]*pIn[idx_f+2];
                
                is3=pIn[idx_f]*(2107.0*pIn[idx_f]-9402.0*pIn[idx_f+1]+7042.0*pIn[idx_f+2]-1854.0*pIn[idx_f+3])
                +pIn[idx_f+1]*(11003.0*pIn[idx_f+1]-17246.0*pIn[idx_f+2]+4642.0*pIn[idx_f+3])
                +pIn[idx_f+2]*(7043.0*pIn[idx_f+2]-3882.0*pIn[idx_f+3])
                +547.0*pIn[idx_f+3]*pIn[idx_f+3];
                
                q40=coeff_weno7_alpha[0][0]*pIn[idx_f-3]+coeff_weno7_alpha[0][1]*pIn[idx_f-2]
                +coeff_weno7_alpha[0][2]*pIn[idx_f-1]+coeff_weno7_alpha[0][3]*pIn[idx_f];
                
                q41=coeff_weno7_alpha[1][0]*pIn[idx_f-2]+coeff_weno7_alpha[1][1]*pIn[idx_f-1]
                +coeff_weno7_alpha[1][2]*pIn[idx_f]+coeff_weno7_alpha[1][3]*pIn[idx_f+1];
                
                q42=coeff_weno7_alpha[2][0]*pIn[idx_f-1]+coeff_weno7_alpha[2][1]*pIn[idx_f]
                +coeff_weno7_alpha[2][2]*pIn[idx_f+1]+coeff_weno7_alpha[2][3]*pIn[idx_f+2];
                
                q43=coeff_weno7_alpha[3][0]*pIn[idx_f]+coeff_weno7_alpha[3][1]*pIn[idx_f+1]
                +coeff_weno7_alpha[3][2]*pIn[idx_f+2]+coeff_weno7_alpha[3][3]*pIn[idx_f+3];
                
                tau7=std::abs(is0-is3);
                alpha0=coeff_weno7_c[0]*(1.0+std::pow(tau7/(ss+is0),2.0));
                alpha1=coeff_weno7_c[1]*(1.0+std::pow(tau7/(ss+is1),2.0));
                alpha2=coeff_weno7_c[2]*(1.0+std::pow(tau7/(ss+is2),2.0));
                alpha3=coeff_weno7_c[3]*(1.0+std::pow(tau7/(ss+is3),2.0));
                
                
                alphaSum=alpha0+alpha1+alpha2+alpha3;
                
                omega0=alpha0/alphaSum;
                omega1=alpha1/alphaSum;
                omega2=alpha2/alphaSum;
                omega3=alpha3/alphaSum;
                
                a[i-ibeg]=0.0;
                b[i-ibeg]=1.0;
                c[i-ibeg]=0.0;
                d[i-ibeg]=omega0*q40+omega1*q41+omega2*q42+omega3*q43;
            }
            
            for(int i=1;i<(iend-ibeg);i++)
            {
                c[i]=c[i]/(b[i]-c[i-1]*a[i]);
                d[i]=(d[i]-d[i-1]*a[i])/(b[i]-c[i-1]*a[i]);
            }
            
            pOut[nx0*ny0*k+nx0*j+iend]=d[iend-ibeg-1];
            
            for(int i=(iend-ibeg-2);i>=0;i--)
            {
                pOut[nx0*ny0*k+nx0*j+i+ibeg+1]=d[i]-c[i]*pOut[nx0*ny0*k+nx0*j+i+ibeg+2];
            }
            
        }
    }
}

void nuc3d::hccs::hccsn(const Field & fieldIN,
                        Field & fieldOUT,
                        const int tilesize)
{
    
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    
    int nx=fieldIN.getSizeX();
    int ny=fieldIN.getSizeY();
    int nz=fieldIN.getSizeZ();
    
    //    vector is very very slow as this function is called many times
    //    std::vector<double> a(nx-2*tilesize+1,0.0);
    //    std::vector<double> b(nx-2*tilesize+1,1.0);
    //    std::vector<double> c(nx-2*tilesize+1,0.0);
    //    std::vector<double> d(nx-2*tilesize+1,0.0);
    double a[256];
    double b[256];
    double c[256];
    double d[256];
    
    int nx0=fieldOUT.getSizeX();
    int ny0=fieldOUT.getSizeY();
    int nz0=fieldOUT.getSizeZ();
    
    int ibeg=(tilesize-1);
    int iend=nx-tilesize;
    int jbeg=0;
    int jend=ny;
    int kbeg=0;
    int kend=nz;
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<(ibeg+1);i++)
            {
                double omega0,omega1,omega2,omega3;
                double alpha0,alpha1,alpha2,alpha3,alphaSum;
                double is0,is1,is2,is3;
                double q40,q41,q42,q43;
                double tau7;
                
                int idx_f=nx*ny*k+nx*j+i;
                int idx_rf=nx0*ny0*k+nx0*j+i+1;
                
                is0=pIn[idx_f+4]*(547.0*pIn[idx_f+4]-3882.0*pIn[idx_f+3]+4642.0*pIn[idx_f+2]-1854.0*pIn[idx_f+1])
                +pIn[idx_f+3]*(7043.0*pIn[idx_f+3]-17246.0*pIn[idx_f+2]+7042.0*pIn[idx_f+1])
                +pIn[idx_f+2]*(11003.0*pIn[idx_f+2]-9402.0*pIn[idx_f+1])
                +2107.0*pIn[idx_f+1]*pIn[idx_f+1];
                
                is1=pIn[idx_f+3]*(267.0*pIn[idx_f+3]-1642.0*pIn[idx_f+2]+1602.0*pIn[idx_f+1]-494.0*pIn[idx_f])
                +pIn[idx_f+2]*(2843.0*pIn[idx_f+2]-5966.0*pIn[idx_f+1]+1922.0*pIn[idx_f])
                +pIn[idx_f+1]*(3443.0*pIn[idx_f+1]-2522.0*pIn[idx_f])
                +547.0*pIn[idx_f]*pIn[idx_f];
                
                
                is2=pIn[idx_f+2]*(547.0*pIn[idx_f+2]-2522.0*pIn[idx_f+1]+1922.0*pIn[idx_f]-494.0*pIn[idx_f-1])
                +pIn[idx_f+1]*(3443.0*pIn[idx_f+1]-5966.0*pIn[idx_f]+1602.0*pIn[idx_f-1])
                +pIn[idx_f]*(2843.0*pIn[idx_f]-1642.0*pIn[idx_f-1])
                +267.0*pIn[idx_f-1]*pIn[idx_f-1];
                
                is3=pIn[idx_f+1]*(2107.0*pIn[idx_f+1]-9402.0*pIn[idx_f]+7042.0*pIn[idx_f-1]-1854.0*pIn[idx_f-2])
                +pIn[idx_f]*(11003.0*pIn[idx_f]-17246.0*pIn[idx_f-1]+4642.0*pIn[idx_f-2])
                +pIn[idx_f-1]*(7043.0*pIn[idx_f-1]-3882.0*pIn[idx_f-2])
                +547.0*pIn[idx_f-2]*pIn[idx_f-2];
                
                q40=coeff_weno7_alpha[0][0]*pIn[idx_f+4]+coeff_weno7_alpha[0][1]*pIn[idx_f+3]
                +coeff_weno7_alpha[0][2]*pIn[idx_f+2]+coeff_weno7_alpha[0][3]*pIn[idx_f+1];
                
                q41=coeff_weno7_alpha[1][0]*pIn[idx_f+3]+coeff_weno7_alpha[1][1]*pIn[idx_f+2]
                +coeff_weno7_alpha[1][2]*pIn[idx_f+1]+coeff_weno7_alpha[1][3]*pIn[idx_f];
                
                q42=coeff_weno7_alpha[2][0]*pIn[idx_f+2]+coeff_weno7_alpha[2][1]*pIn[idx_f+1]
                +coeff_weno7_alpha[2][2]*pIn[idx_f]+coeff_weno7_alpha[2][3]*pIn[idx_f-1];
                
                q43=coeff_weno7_alpha[3][0]*pIn[idx_f+1]+coeff_weno7_alpha[3][1]*pIn[idx_f]
                +coeff_weno7_alpha[3][2]*pIn[idx_f-1]+coeff_weno7_alpha[3][3]*pIn[idx_f-2];
                
                tau7=std::abs(is0-is3);
                alpha0=coeff_weno7_c[0]*(1.0+std::pow(tau7/(ss+is0),2.0));
                alpha1=coeff_weno7_c[1]*(1.0+std::pow(tau7/(ss+is1),2.0));
                alpha2=coeff_weno7_c[2]*(1.0+std::pow(tau7/(ss+is2),2.0));
                alpha3=coeff_weno7_c[3]*(1.0+std::pow(tau7/(ss+is3),2.0));
                
                
                alphaSum=alpha0+alpha1+alpha2+alpha3;
                
                omega0=alpha0/alphaSum;
                omega1=alpha1/alphaSum;
                omega2=alpha2/alphaSum;
                omega3=alpha3/alphaSum;
                
                a[i-ibeg]=0.0;
                b[i-ibeg]=1.0;
                c[i-ibeg]=0.0;
                d[i-ibeg]=omega0*q40+omega1*q41+omega2*q42+omega3*q43;
            }
            
            for(int i=ibeg+1;i<(iend-1);i++)
            {
                double omega0,omega1,omega2;
                double alpha0,alpha1,alpha2,alphaSum;
                double is0,is1,is2;
                double q7;
                
                int idx_f=nx*ny*k+nx*j+i;
                
                is0= coeff_weno5_gamma0*std::pow((pIn[idx_f+3]-2.0*pIn[idx_f+2]+pIn[idx_f+1]),2)
                +coeff_weno5_gamma1*std::pow((pIn[idx_f+3]-4.0*pIn[idx_f+2]+3.0*pIn[idx_f+1]),2);
                
                is1= coeff_weno5_gamma0*std::pow((pIn[idx_f+2]-2.0*pIn[idx_f+1]+pIn[idx_f]),2)
                +coeff_weno5_gamma1*std::pow((pIn[idx_f+2]-pIn[idx_f]),2);
                
                is2= coeff_weno5_gamma0*std::pow((pIn[idx_f+1]-2.0*pIn[idx_f]+pIn[idx_f-1]),2)
                +coeff_weno5_gamma1*std::pow((3.0*pIn[idx_f+1]-4.0*pIn[idx_f]+pIn[idx_f-1]),2);
                
                q7=(-pIn[idx_f+3]+19.0*pIn[idx_f+2]+239.0*pIn[idx_f+1]+159.0*pIn[idx_f]+4.0*pIn[idx_f-1])/420.0;
                
                double tau=std::abs(is2-is0);
                
                alpha0=coeff_crweno5_c[0]*(1.0+std::pow(tau/(is0+ss),p));
                alpha1=coeff_crweno5_c[1]*(1.0+std::pow(tau/(is1+ss),p));
                alpha2=coeff_crweno5_c[2]*(1.0+std::pow(tau/(is2+ss),p));
                
                alphaSum=alpha0+alpha1+alpha2;
                
                omega0=alpha0/alphaSum;
                omega1=alpha1/alphaSum;
                omega2=alpha2/alphaSum;
                
                double theta=1.0/(1.0+std::pow(alphaSum-1.0,2.0));
                
                c[i-ibeg]=(2.0/3.0*omega0+1.0/3.0*omega1)*(1.0-theta)
                +2.0/7.0*theta;
                b[i-ibeg]=(1.0/3.0*omega0+2.0/3.0*(omega1+omega2))*(1.0-theta)
                +4.0/7.0*theta;
                a[i-ibeg]=(1.0/3.0*omega2)*(1.0-theta)
                +1.0/7.0*theta;
                d[i-ibeg]=(1.0/6.0*omega0*pIn[idx_f+2]
                           +(5.0*(omega0+omega1)+omega2)*pIn[idx_f+1]/6.0
                           +(omega1+5.0*omega2)*pIn[idx_f]/6.0)*(1.0-theta)+q7*theta;
            }
            
            for(int i=iend-1;i<iend;i++)
            {
                double omega0,omega1,omega2,omega3;
                double alpha0,alpha1,alpha2,alpha3,alphaSum;
                double is0,is1,is2,is3;
                double q40,q41,q42,q43;
                double tau7;
                
                int idx_f=nx*ny*k+nx*j+i;
                
                is0=pIn[idx_f+4]*(547.0*pIn[idx_f+4]-3882.0*pIn[idx_f+3]+4642.0*pIn[idx_f+2]-1854.0*pIn[idx_f+1])
                +pIn[idx_f+3]*(7043.0*pIn[idx_f+3]-17246.0*pIn[idx_f+2]+7042.0*pIn[idx_f+1])
                +pIn[idx_f+2]*(11003.0*pIn[idx_f+2]-9402.0*pIn[idx_f+1])
                +2107.0*pIn[idx_f+1]*pIn[idx_f+1];
                
                is1=pIn[idx_f+3]*(267.0*pIn[idx_f+3]-1642.0*pIn[idx_f+2]+1602.0*pIn[idx_f+1]-494.0*pIn[idx_f])
                +pIn[idx_f+2]*(2843.0*pIn[idx_f+2]-5966.0*pIn[idx_f+1]+1922.0*pIn[idx_f])
                +pIn[idx_f+1]*(3443.0*pIn[idx_f+1]-2522.0*pIn[idx_f])
                +547.0*pIn[idx_f]*pIn[idx_f];
                
                
                is2=pIn[idx_f+2]*(547.0*pIn[idx_f+2]-2522.0*pIn[idx_f+1]+1922.0*pIn[idx_f]-494.0*pIn[idx_f-1])
                +pIn[idx_f+1]*(3443.0*pIn[idx_f+1]-5966.0*pIn[idx_f]+1602.0*pIn[idx_f-1])
                +pIn[idx_f]*(2843.0*pIn[idx_f]-1642.0*pIn[idx_f-1])
                +267.0*pIn[idx_f-1]*pIn[idx_f-1];
                
                is3=pIn[idx_f+1]*(2107.0*pIn[idx_f+1]-9402.0*pIn[idx_f]+7042.0*pIn[idx_f-1]-1854.0*pIn[idx_f-2])
                +pIn[idx_f]*(11003.0*pIn[idx_f]-17246.0*pIn[idx_f-1]+4642.0*pIn[idx_f-2])
                +pIn[idx_f-1]*(7043.0*pIn[idx_f-1]-3882.0*pIn[idx_f-2])
                +547.0*pIn[idx_f-2]*pIn[idx_f-2];
                
                q40=coeff_weno7_alpha[0][0]*pIn[idx_f+4]+coeff_weno7_alpha[0][1]*pIn[idx_f+3]
                +coeff_weno7_alpha[0][2]*pIn[idx_f+2]+coeff_weno7_alpha[0][3]*pIn[idx_f+1];
                
                q41=coeff_weno7_alpha[1][0]*pIn[idx_f+3]+coeff_weno7_alpha[1][1]*pIn[idx_f+2]
                +coeff_weno7_alpha[1][2]*pIn[idx_f+1]+coeff_weno7_alpha[1][3]*pIn[idx_f];
                
                q42=coeff_weno7_alpha[2][0]*pIn[idx_f+2]+coeff_weno7_alpha[2][1]*pIn[idx_f+1]
                +coeff_weno7_alpha[2][2]*pIn[idx_f]+coeff_weno7_alpha[2][3]*pIn[idx_f-1];
                
                q43=coeff_weno7_alpha[3][0]*pIn[idx_f+1]+coeff_weno7_alpha[3][1]*pIn[idx_f]
                +coeff_weno7_alpha[3][2]*pIn[idx_f-1]+coeff_weno7_alpha[3][3]*pIn[idx_f-2];
                
                tau7=std::abs(is0-is3);
                alpha0=coeff_weno7_c[0]*(1.0+std::pow(tau7/(ss+is0),2.0));
                alpha1=coeff_weno7_c[1]*(1.0+std::pow(tau7/(ss+is1),2.0));
                alpha2=coeff_weno7_c[2]*(1.0+std::pow(tau7/(ss+is2),2.0));
                alpha3=coeff_weno7_c[3]*(1.0+std::pow(tau7/(ss+is3),2.0));
                
                
                alphaSum=alpha0+alpha1+alpha2+alpha3;
                
                omega0=alpha0/alphaSum;
                omega1=alpha1/alphaSum;
                omega2=alpha2/alphaSum;
                omega3=alpha3/alphaSum;
                
                a[i-ibeg]=0.0;
                b[i-ibeg]=1.0;
                c[i-ibeg]=0.0;
                d[i-ibeg]=omega0*q40+omega1*q41+omega2*q42+omega3*q43;
            }
            
            for(int i=1;i<(iend-ibeg);i++)
            {
                c[i]=c[i]/(b[i]-c[i-1]*a[i]);
                d[i]=(d[i]-d[i-1]*a[i])/(b[i]-c[i-1]*a[i]);
            }
            
            pOut[nx0*ny0*k+nx0*j+iend]=d[iend-ibeg-1];
            
            for(int i=(iend-ibeg-2);i>=0;i--)
            {
                pOut[nx0*ny0*k+nx0*j+i+ibeg+1]=d[i]-c[i]*pOut[nx0*ny0*k+nx0*j+i+ibeg+2];
            }
        }
    }
}

void nuc3d::hccs::interpolationBoundaryL(const Field & fieldIN,
                                         const Field & boundaryL,
                                         const int uw,//uw=(1,-1)
                                         Field & fieldOUT,
                                         const int tilesize)
{
    switch(uw)
    {
        case 1:
            hccspBL(fieldIN,boundaryL, fieldOUT, tilesize);
            break;
        case -1:
            hccsnBL(fieldIN,boundaryL, fieldOUT, tilesize);
            break;
        default:
            std::cout<<"weno error: no such direction"<<std::endl;
            exit(-1);
    }
    
}

void nuc3d::hccs::hccspBL(const Field & fieldIN,
                          const Field & boundaryL,
                          Field & fieldOUT,
                          const int tilesize)
{
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    double *pBND=boundaryL.getDataPtr();
    
    double omega0,omega1,omega2,omega3;
    double alpha0,alpha1,alpha2,alpha3,alphaSum;
    double is0,is1,is2,is3;
    double q40,q41,q42,q43;
    double tau7;
    
    int nx=fieldIN.getSizeX();
    int ny=fieldIN.getSizeY();
    int nz=fieldIN.getSizeZ();
    
    int nx0=fieldOUT.getSizeX();
    int ny0=fieldOUT.getSizeY();
    int nz0=fieldOUT.getSizeZ();
    
    int nxBND=boundaryL.getSizeX();
    int nyBND=boundaryL.getSizeY();
    int nzBND=boundaryL.getSizeZ();
    
    int ibeg=0;
    int iend=tilesize;
    int jbeg=0;
    int jend=ny;
    int kbeg=0;
    int kend=nz;
    
    double f[7];
    
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i;
                for(int z=-3;z<=3;z++)
                {
                    if((i+z-1)<0)
                    {
                        int idx_BND=nxBND*nyBND*k+nxBND*j+(i+z-1+nxBND);
                        
                        f[3+z]=pBND[idx_BND];
                    }
                    else
                    {
                        int idx_f=nx*ny*k+nx*j+i+z-1;
                        f[3+z]=pIn[idx_f];
                    }
                }
                
                
                is0=f[0]*(547.0*f[0]-3882.0*f[1]+4642.0*f[2]-1854.0*f[3])
                +f[1]*(7043.0*f[1]-17246.0*f[2]+7042.0*f[3])
                +f[2]*(11003.0*f[2]-9402.0*f[3])
                +2107.0*f[3]*f[3];
                
                is1=f[1]*(267.0*f[1]-1642.0*f[2]+1602.0*f[3]-494.0*f[4])
                +f[2]*(2843.0*f[2]-5966.0*f[3]+1922.0*f[4])
                +f[3]*(3443.0*f[3]-2522.0*f[4])
                +547.0*f[4]*f[4];
                
                
                is2=f[2]*(547.0*f[2]-2522.0*f[3]+1922.0*f[4]-494.0*f[5])
                +f[3]*(3443.0*f[3]-5966.0*f[4]+1602.0*f[5])
                +f[4]*(2843.0*f[4]-1642.0*f[5])
                +267.0*f[5]*f[5];
                
                is3=f[3]*(2107.0*f[3]-9402.0*f[4]+7042.0*f[5]-1854.0*f[6])
                +f[4]*(11003.0*f[4]-17246.0*f[5]+4642.0*f[6])
                +f[5]*(7043.0*f[5]-3882.0*f[6])
                +547.0*f[6]*f[6];
                
                q40=coeff_weno7_alpha[0][0]*f[0]+coeff_weno7_alpha[0][1]*f[1]
                +coeff_weno7_alpha[0][2]*f[2]+coeff_weno7_alpha[0][3]*f[3];
                
                q41=coeff_weno7_alpha[1][0]*f[1]+coeff_weno7_alpha[1][1]*f[2]
                +coeff_weno7_alpha[1][2]*f[3]+coeff_weno7_alpha[1][3]*f[4];
                
                q42=coeff_weno7_alpha[2][0]*f[2]+coeff_weno7_alpha[2][1]*f[3]
                +coeff_weno7_alpha[2][2]*f[4]+coeff_weno7_alpha[2][3]*f[5];
                
                q43=coeff_weno7_alpha[3][0]*f[3]+coeff_weno7_alpha[3][1]*f[4]
                +coeff_weno7_alpha[3][2]*f[5]+coeff_weno7_alpha[3][3]*f[6];
                
                tau7=std::abs(is0-is3);
                alpha0=coeff_weno7_c[0]*(1.0+std::pow(tau7/(ss+is0),2.0));
                alpha1=coeff_weno7_c[1]*(1.0+std::pow(tau7/(ss+is1),2.0));
                alpha2=coeff_weno7_c[2]*(1.0+std::pow(tau7/(ss+is2),2.0));
                alpha3=coeff_weno7_c[3]*(1.0+std::pow(tau7/(ss+is3),2.0));
                
                
                alphaSum=alpha0+alpha1+alpha2+alpha3;
                
                omega0=alpha0/alphaSum;
                omega1=alpha1/alphaSum;
                omega2=alpha2/alphaSum;
                omega3=alpha3/alphaSum;
                
                pOut[idx_rf]=omega0*q40+omega1*q41+omega2*q42+omega3*q43;
            }
        }
    }
}

void nuc3d::hccs::hccsnBL(const Field & fieldIN,
                          const Field & boundaryL,
                          Field & fieldOUT,
                          const int tilesize)
{
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    double *pBND=boundaryL.getDataPtr();
    
    double omega0,omega1,omega2,omega3;
    double alpha0,alpha1,alpha2,alpha3,alphaSum;
    double is0,is1,is2,is3;
    double q40,q41,q42,q43;
    double tau7;
    
    int nx=fieldIN.getSizeX();
    int ny=fieldIN.getSizeY();
    int nz=fieldIN.getSizeZ();
    
    int nx0=fieldOUT.getSizeX();
    int ny0=fieldOUT.getSizeY();
    int nz0=fieldOUT.getSizeZ();
    
    int nxBND=boundaryL.getSizeX();
    int nyBND=boundaryL.getSizeY();
    int nzBND=boundaryL.getSizeZ();
    
    int ibeg=0;
    int iend=tilesize;
    int jbeg=0;
    int jend=ny;
    int kbeg=0;
    int kend=nz;
    
    double f[7];
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i;
                
                for(int z=-3;z<=3;z++)
                {
                    if((i-z)<0)
                    {
                        int idx_BND=nxBND*nyBND*k+nxBND*j+(i-z+nxBND);
                        
                        f[3+z]=pBND[idx_BND];
                    }
                    else
                    {
                        int idx_f=nx*ny*k+nx*j+i-z;
                        f[3+z]=pIn[idx_f];
                    }
                }
                
                
                
                is0=f[0]*(547.0*f[0]-3882.0*f[1]+4642.0*f[2]-1854.0*f[3])
                +f[1]*(7043.0*f[1]-17246.0*f[2]+7042.0*f[3])
                +f[2]*(11003.0*f[2]-9402.0*f[3])
                +2107.0*f[3]*f[3];
                
                is1=f[1]*(267.0*f[1]-1642.0*f[2]+1602.0*f[3]-494.0*f[4])
                +f[2]*(2843.0*f[2]-5966.0*f[3]+1922.0*f[4])
                +f[3]*(3443.0*f[3]-2522.0*f[4])
                +547.0*f[4]*f[4];
                
                
                is2=f[2]*(547.0*f[2]-2522.0*f[3]+1922.0*f[4]-494.0*f[5])
                +f[3]*(3443.0*f[3]-5966.0*f[4]+1602.0*f[5])
                +f[4]*(2843.0*f[4]-1642.0*f[5])
                +267.0*f[5]*f[5];
                
                is3=f[3]*(2107.0*f[3]-9402.0*f[4]+7042.0*f[5]-1854.0*f[6])
                +f[4]*(11003.0*f[4]-17246.0*f[5]+4642.0*f[6])
                +f[5]*(7043.0*f[5]-3882.0*f[6])
                +547.0*f[6]*f[6];
                
                q40=coeff_weno7_alpha[0][0]*f[0]+coeff_weno7_alpha[0][1]*f[1]
                +coeff_weno7_alpha[0][2]*f[2]+coeff_weno7_alpha[0][3]*f[3];
                
                q41=coeff_weno7_alpha[1][0]*f[1]+coeff_weno7_alpha[1][1]*f[2]
                +coeff_weno7_alpha[1][2]*f[3]+coeff_weno7_alpha[1][3]*f[4];
                
                q42=coeff_weno7_alpha[2][0]*f[2]+coeff_weno7_alpha[2][1]*f[3]
                +coeff_weno7_alpha[2][2]*f[4]+coeff_weno7_alpha[2][3]*f[5];
                
                q43=coeff_weno7_alpha[3][0]*f[3]+coeff_weno7_alpha[3][1]*f[4]
                +coeff_weno7_alpha[3][2]*f[5]+coeff_weno7_alpha[3][3]*f[6];
                
                tau7=std::abs(is0-is3);
                alpha0=coeff_weno7_c[0]*(1.0+std::pow(tau7/(ss+is0),2.0));
                alpha1=coeff_weno7_c[1]*(1.0+std::pow(tau7/(ss+is1),2.0));
                alpha2=coeff_weno7_c[2]*(1.0+std::pow(tau7/(ss+is2),2.0));
                alpha3=coeff_weno7_c[3]*(1.0+std::pow(tau7/(ss+is3),2.0));
                
                
                alphaSum=alpha0+alpha1+alpha2+alpha3;
                
                omega0=alpha0/alphaSum;
                omega1=alpha1/alphaSum;
                omega2=alpha2/alphaSum;
                omega3=alpha3/alphaSum;
                
                pOut[idx_rf]=omega0*q40+omega1*q41+omega2*q42+omega3*q43;
            }
        }
    }
}

void nuc3d::hccs::interpolationBoundaryR(const Field & fieldIN,
                                         const Field & boundaryR,
                                         const int uw,//uw=(1,-1)
                                         Field & fieldOUT,
                                         const int tilesize)
{
    switch(uw)
    {
        case 1:
            hccspBR(fieldIN,boundaryR, fieldOUT, tilesize);
            break;
        case -1:
            hccsnBR(fieldIN,boundaryR, fieldOUT, tilesize);
            break;
        default:
            std::cout<<"weno error: no such direction"<<std::endl;
            exit(-1);
    }
    
}


void nuc3d::hccs::hccspBR(const Field & fieldIN,
                          const Field & boundaryR,
                          Field & fieldOUT,
                          const int tilesize)
{
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    double *pBND=boundaryR.getDataPtr();
    
    double omega0,omega1,omega2,omega3;
    double alpha0,alpha1,alpha2,alpha3,alphaSum;
    double is0,is1,is2,is3;
    double q40,q41,q42,q43;
    double tau7;
    
    int nx=fieldIN.getSizeX();
    int ny=fieldIN.getSizeY();
    int nz=fieldIN.getSizeZ();
    
    int nx0=fieldOUT.getSizeX();
    int ny0=fieldOUT.getSizeY();
    int nz0=fieldOUT.getSizeZ();
    
    int nxBND=boundaryR.getSizeX();
    int nyBND=boundaryR.getSizeY();
    int nzBND=boundaryR.getSizeZ();
    
    int ibeg=nx-tilesize;
    int iend=nx;
    int jbeg=0;
    int jend=ny;
    int kbeg=0;
    int kend=nz;
    
    double f[7];
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i+1;
                for(int z=-3;z<=3;z++)
                {
                    if((i+z)>=nx)
                    {
                        int idx_BND=nxBND*nyBND*k+nxBND*j+i+z-nx;
                        
                        f[3+z]=pBND[idx_BND];
                    }
                    else
                    {
                        int idx_f=nx*ny*k+nx*j+i+z;
                        f[3+z]=pIn[idx_f];
                    }
                }
                
                
                
                is0=f[0]*(547.0*f[0]-3882.0*f[1]+4642.0*f[2]-1854.0*f[3])
                +f[1]*(7043.0*f[1]-17246.0*f[2]+7042.0*f[3])
                +f[2]*(11003.0*f[2]-9402.0*f[3])
                +2107.0*f[3]*f[3];
                
                is1=f[1]*(267.0*f[1]-1642.0*f[2]+1602.0*f[3]-494.0*f[4])
                +f[2]*(2843.0*f[2]-5966.0*f[3]+1922.0*f[4])
                +f[3]*(3443.0*f[3]-2522.0*f[4])
                +547.0*f[4]*f[4];
                
                
                is2=f[2]*(547.0*f[2]-2522.0*f[3]+1922.0*f[4]-494.0*f[5])
                +f[3]*(3443.0*f[3]-5966.0*f[4]+1602.0*f[5])
                +f[4]*(2843.0*f[4]-1642.0*f[5])
                +267.0*f[5]*f[5];
                
                is3=f[3]*(2107.0*f[3]-9402.0*f[4]+7042.0*f[5]-1854.0*f[6])
                +f[4]*(11003.0*f[4]-17246.0*f[5]+4642.0*f[6])
                +f[5]*(7043.0*f[5]-3882.0*f[6])
                +547.0*f[6]*f[6];
                
                q40=coeff_weno7_alpha[0][0]*f[0]+coeff_weno7_alpha[0][1]*f[1]
                +coeff_weno7_alpha[0][2]*f[2]+coeff_weno7_alpha[0][3]*f[3];
                
                q41=coeff_weno7_alpha[1][0]*f[1]+coeff_weno7_alpha[1][1]*f[2]
                +coeff_weno7_alpha[1][2]*f[3]+coeff_weno7_alpha[1][3]*f[4];
                
                q42=coeff_weno7_alpha[2][0]*f[2]+coeff_weno7_alpha[2][1]*f[3]
                +coeff_weno7_alpha[2][2]*f[4]+coeff_weno7_alpha[2][3]*f[5];
                
                q43=coeff_weno7_alpha[3][0]*f[3]+coeff_weno7_alpha[3][1]*f[4]
                +coeff_weno7_alpha[3][2]*f[5]+coeff_weno7_alpha[3][3]*f[6];
                
                tau7=std::abs(is0-is3);
                alpha0=coeff_weno7_c[0]*(1.0+std::pow(tau7/(ss+is0),2.0));
                alpha1=coeff_weno7_c[1]*(1.0+std::pow(tau7/(ss+is1),2.0));
                alpha2=coeff_weno7_c[2]*(1.0+std::pow(tau7/(ss+is2),2.0));
                alpha3=coeff_weno7_c[3]*(1.0+std::pow(tau7/(ss+is3),2.0));
                
                
                alphaSum=alpha0+alpha1+alpha2+alpha3;
                
                omega0=alpha0/alphaSum;
                omega1=alpha1/alphaSum;
                omega2=alpha2/alphaSum;
                omega3=alpha3/alphaSum;
                
                pOut[idx_rf]=omega0*q40+omega1*q41+omega2*q42+omega3*q43;
            }
        }
    }
}

void nuc3d::hccs::hccsnBR(const Field & fieldIN,
                          const Field & boundaryR,
                          Field & fieldOUT,
                          const int tilesize)
{
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    double *pBND=boundaryR.getDataPtr();
    
    double omega0,omega1,omega2,omega3;
    double alpha0,alpha1,alpha2,alpha3,alphaSum;
    double is0,is1,is2,is3;
    double q40,q41,q42,q43;
    double tau7;
    
    int nx=fieldIN.getSizeX();
    int ny=fieldIN.getSizeY();
    int nz=fieldIN.getSizeZ();
    
    int nx0=fieldOUT.getSizeX();
    int ny0=fieldOUT.getSizeY();
    int nz0=fieldOUT.getSizeZ();
    
    int nxBND=boundaryR.getSizeX();
    int nyBND=boundaryR.getSizeY();
    int nzBND=boundaryR.getSizeZ();
    
    int ibeg=nx-tilesize;
    int iend=nx;
    int jbeg=0;
    int jend=ny;
    int kbeg=0;
    int kend=nz;
    
    double f[7];
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i+1;
                
                for(int z=-3;z<=3;z++)
                {
                    if((i-z+1)>=nx)
                    {
                        int idx_BND=nxBND*nyBND*k+nxBND*j+i-z+1-nx;
                        
                        f[3+z]=pBND[idx_BND];
                    }
                    else
                    {
                        int idx_f=nx*ny*k+nx*j+i-z+1;
                        f[3+z]=pIn[idx_f];
                    }
                }
                
                
                
                is0=f[0]*(547.0*f[0]-3882.0*f[1]+4642.0*f[2]-1854.0*f[3])
                +f[1]*(7043.0*f[1]-17246.0*f[2]+7042.0*f[3])
                +f[2]*(11003.0*f[2]-9402.0*f[3])
                +2107.0*f[3]*f[3];
                
                is1=f[1]*(267.0*f[1]-1642.0*f[2]+1602.0*f[3]-494.0*f[4])
                +f[2]*(2843.0*f[2]-5966.0*f[3]+1922.0*f[4])
                +f[3]*(3443.0*f[3]-2522.0*f[4])
                +547.0*f[4]*f[4];
                
                
                is2=f[2]*(547.0*f[2]-2522.0*f[3]+1922.0*f[4]-494.0*f[5])
                +f[3]*(3443.0*f[3]-5966.0*f[4]+1602.0*f[5])
                +f[4]*(2843.0*f[4]-1642.0*f[5])
                +267.0*f[5]*f[5];
                
                is3=f[3]*(2107.0*f[3]-9402.0*f[4]+7042.0*f[5]-1854.0*f[6])
                +f[4]*(11003.0*f[4]-17246.0*f[5]+4642.0*f[6])
                +f[5]*(7043.0*f[5]-3882.0*f[6])
                +547.0*f[6]*f[6];
                
                q40=coeff_weno7_alpha[0][0]*f[0]+coeff_weno7_alpha[0][1]*f[1]
                +coeff_weno7_alpha[0][2]*f[2]+coeff_weno7_alpha[0][3]*f[3];
                
                q41=coeff_weno7_alpha[1][0]*f[1]+coeff_weno7_alpha[1][1]*f[2]
                +coeff_weno7_alpha[1][2]*f[3]+coeff_weno7_alpha[1][3]*f[4];
                
                q42=coeff_weno7_alpha[2][0]*f[2]+coeff_weno7_alpha[2][1]*f[3]
                +coeff_weno7_alpha[2][2]*f[4]+coeff_weno7_alpha[2][3]*f[5];
                
                q43=coeff_weno7_alpha[3][0]*f[3]+coeff_weno7_alpha[3][1]*f[4]
                +coeff_weno7_alpha[3][2]*f[5]+coeff_weno7_alpha[3][3]*f[6];
                
                tau7=std::abs(is0-is3);
                alpha0=coeff_weno7_c[0]*(1.0+std::pow(tau7/(ss+is0),2.0));
                alpha1=coeff_weno7_c[1]*(1.0+std::pow(tau7/(ss+is1),2.0));
                alpha2=coeff_weno7_c[2]*(1.0+std::pow(tau7/(ss+is2),2.0));
                alpha3=coeff_weno7_c[3]*(1.0+std::pow(tau7/(ss+is3),2.0));
                
                
                alphaSum=alpha0+alpha1+alpha2+alpha3;
                
                omega0=alpha0/alphaSum;
                omega1=alpha1/alphaSum;
                omega2=alpha2/alphaSum;
                omega3=alpha3/alphaSum;
                
                pOut[idx_rf]=omega0*q40+omega1*q41+omega2*q42+omega3*q43;
            }
        }
    }
}
