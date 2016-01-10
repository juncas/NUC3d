//
//  reconstructionBoundaryScheme.cpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/10.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//

#include "reconstructionBoundaryScheme.hpp"

#include "schemes.hpp"

nuc3d::ReconstructionboundaryScheme::ReconstructionboundaryScheme()
{
    
}

nuc3d::ReconstructionboundaryScheme::~ReconstructionboundaryScheme()
{
    
}

void nuc3d::ReconstructionboundaryScheme::interpolationInner(const Field & fieldIN,
                                                             const int uw,//uw=(1,-1)
                                                             Field & fieldOUT,
                                                             const int tilesize)
{
    
}

void nuc3d::ReconstructionboundaryScheme::interpolationBoundaryL(const Field & fieldIN,
                                                                 const Field & boundaryL,
                                                                 const int uw,//uw=(1,-1)
                                                                 Field & fieldOUT,
                                                                 const int tilesize)
{
    switch(uw)
    {
        case 1:
            boundary_pL(fieldIN,boundaryL, fieldOUT, tilesize);
            break;
        case -1:
            boundary_nL(fieldIN,boundaryL, fieldOUT, tilesize);
            break;
        default:
            std::cout<<"weno error: no such direction"<<std::endl;
            exit(-1);
    }
    
}

void nuc3d::ReconstructionboundaryScheme::boundary_pL(const Field & fieldIN,
                                                      const Field & boundaryL,
                                                      Field & fieldOUT,
                                                      const int tilesize)
{
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    double *pBND=boundaryL.getDataPtr();
    
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
    
    double *f1st;
    double *f2nd;
    double *f3rd;
    double *f5th;
    
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            int idx_rf=nx0*ny0*k+nx0*j;
            int idx_f=nx*ny*k+nx*j;
            int idx_BND=nxBND*nyBND*k+nxBND*j+nxBND-1;
            
            f1st=pBND+idx_BND;
            f2nd=pIn+idx_f;
            f3rd=pIn+idx_f;
            f5th=pIn+idx_f;
            
            firstOrderP(f1st,pOut+idx_rf);
            secondOrderP(f2nd,pOut+idx_rf+1);
            thirdOrderP(f3rd,pOut+idx_rf+2);
            fifthOrderP(f5th,pOut+idx_rf+3);
            
        }
    }
}

void nuc3d::ReconstructionboundaryScheme::boundary_nL(const Field & fieldIN,
                                                      const Field & boundaryL,
                                                      Field & fieldOUT,
                                                      const int tilesize)
{
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    double *pBND=boundaryL.getDataPtr();
    
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
    
    double *f1st;
    double *f2nd;
    double *f3rd;
    double *f5th;
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            int idx_rf=nx0*ny0*k+nx0*j;
            int idx_f=nx*ny*k+nx*j;
            
            f1st=pIn+idx_f;
            f2nd=pIn+idx_f+1;
            f3rd=pIn+idx_f+2;
            f5th=pIn+idx_f+3;
            
            firstOrderP(f1st,pOut+idx_rf);
            secondOrderP(f2nd,pOut+idx_rf+1);
            thirdOrderP(f3rd,pOut+idx_rf+2);
            fifthOrderP(f5th,pOut+idx_rf+3);
        }
    }
}

void nuc3d::ReconstructionboundaryScheme::interpolationBoundaryR(const Field & fieldIN,
                                                                 const Field & boundaryR,
                                                                 const int uw,//uw=(1,-1)
                                                                 Field & fieldOUT,
                                                                 const int tilesize)
{
    switch(uw)
    {
        case 1:
            boundary_pR(fieldIN,boundaryR, fieldOUT, tilesize);
            break;
        case -1:
            boundary_nR(fieldIN,boundaryR, fieldOUT, tilesize);
            break;
        default:
            std::cout<<"weno error: no such direction"<<std::endl;
            exit(-1);
    }
    
}


void nuc3d::ReconstructionboundaryScheme::boundary_pR(const Field & fieldIN,
                                                      const Field & boundaryR,
                                                      Field & fieldOUT,
                                                      const int tilesize)
{
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    double *pBND=boundaryR.getDataPtr();
    
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
    
    double *f1st;
    double *f2nd;
    double *f3rd;
    double *f5th;
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            int idx_rf=nx0*ny0*k+nx0*j+nx0-1;
            int idx_f=nx*ny*k+nx*j+nx-1;
            
            f1st=pIn+idx_f;
            f2nd=pIn+idx_f-1;
            f3rd=pIn+idx_f-3;
            f5th=pIn+idx_f-5;
            
            firstOrderP(f1st,pOut+idx_rf);
            secondOrderP(f2nd,pOut+idx_rf-1);
            thirdOrderP(f3rd,pOut+idx_rf-2);
            fifthOrderP(f5th,pOut+idx_rf-3);
        }
    }
}

void nuc3d::ReconstructionboundaryScheme::boundary_nR(const Field & fieldIN,
                                                      const Field & boundaryR,
                                                      Field & fieldOUT,
                                                      const int tilesize)
{
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    double *pBND=boundaryR.getDataPtr();
    
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
    
    double *f1st;
    double *f2nd;
    double *f3rd;
    double *f5th;
    
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            int idx_rf=nx0*ny0*k+nx0*j+nx0-1;
            int idx_f=nx*ny*k+nx*j+nx-1;
            int idx_BND=nxBND*nyBND*k+nxBND*j;
            
            f1st=pBND+idx_BND;
            f2nd=pIn+idx_f;
            f3rd=pIn+idx_f-1;
            f5th=pIn+idx_f-2;
            
            firstOrderP(f1st,pOut+idx_rf);
            secondOrderP(f2nd,pOut+idx_rf-1);
            thirdOrderP(f3rd,pOut+idx_rf-2);
            fifthOrderP(f5th,pOut+idx_rf-3);
        }
    }
}


void nuc3d::ReconstructionboundaryScheme::firstOrderP(double *f,double *fout)
{
    
    *fout=*f;
    
}

void nuc3d::ReconstructionboundaryScheme::secondOrderP(double *f,double *fout)
{
    *fout=0.5*(f[0]+f[1]);
}

void nuc3d::ReconstructionboundaryScheme::thirdOrderP(double *f,double *fout)
{
    double is0=pow((f[0]-f[1]),2);
    double is1=pow((f[1]-f[2]),2);
    
    
    double q20=-0.5*f[0]+1.5*f[1];
    double q21=0.5*f[1]+0.5*f[2];
    
    double aa0=1.0/3.0/pow(1.0e-6+is0,2);
    double aa1=2.0/3.0/pow(1.0e-6+is1,2);
    
    double w0=aa0/(aa0+aa1);
    double w1=aa1/(aa0+aa1);
    
    *fout=w0*q20+w1*q21;
}

void nuc3d::ReconstructionboundaryScheme::fifthOrderP(double *f,double *fout)
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
    
    alpha0=coeff_weno5_c[0]/pow((1.0e-6+is0),2);
    alpha1=coeff_weno5_c[1]/pow((1.0e-6+is1),2);
    alpha2=coeff_weno5_c[2]/pow((1.0e-6+is2),2);
    
    
    alphaSum=alpha0+alpha1+alpha2;
    
    omega0=alpha0/alphaSum;
    omega1=alpha1/alphaSum;
    omega2=alpha2/alphaSum;
    
    *fout=omega0*q30+omega1*q31+omega2*q32;
    
}

void nuc3d::ReconstructionboundaryScheme::firstOrderN(double *f,double *fout)
{
    *fout=*f;
}

void nuc3d::ReconstructionboundaryScheme::secondOrderN(double *f,double *fout)
{
    double f0=*f;
    double f1=*(f-1);
    
    *fout=0.5*(f0+f1);
}

void nuc3d::ReconstructionboundaryScheme::thirdOrderN(double *f,double *fout)
{
    double f0=*(f+1);
    double f1=*f;
    double f2=*(f-1);
    
    double is0=pow((f0-f1),2);
    double is1=pow((f1-f2),2);
    
    
    double q20=-0.5*f0+1.5*f1;
    double q21=0.5*f1+0.5*f2;
    
    double aa0=1.0/3.0/pow(1.0e-6+is0,2);
    double aa1=2.0/3.0/pow(1.0e-6+is1,2);
    
    double w0=aa0/(aa0+aa1);
    double w1=aa1/(aa0+aa1);
    
    *fout=w0*q20+w1*q21;
}

void nuc3d::ReconstructionboundaryScheme::fifthOrderN(double *f,double *fout)
{
    
    double omega0,omega1,omega2;
    double alpha0,alpha1,alpha2,alphaSum;
    double is0,is1,is2;
    double q30,q31,q32;
    
    double f0=*(f+2);
    double f1=*(f+1);
    double f2=*f;
    double f3=*(f-1);
    double f4=*(f-2);
    
    is0= coeff_weno5_gamma0*pow((    f0-2.0*f1+    f2),2)
    +coeff_weno5_gamma1*pow((    f0-4.0*f1+3.0*f2),2);
    
    is1= coeff_weno5_gamma0*pow((    f1-2.0*f2+    f3),2)
    +coeff_weno5_gamma1*pow((    f1-             f3),2);
    
    is2= coeff_weno5_gamma0*pow((    f2-2.0*f3+    f4),2)
    +coeff_weno5_gamma1*pow((3.0*f2-4.0*f3+    f4),2);
    
    
    q30= coeff_weno5_alpha[0][0]*f0
    +coeff_weno5_alpha[0][1]*f1
    +coeff_weno5_alpha[0][2]*f2;
    
    q31= coeff_weno5_alpha[1][0]*f1
    +coeff_weno5_alpha[1][1]*f2
    +coeff_weno5_alpha[1][2]*f3;
    
    q32= coeff_weno5_alpha[2][0]*f2
    +coeff_weno5_alpha[2][1]*f3
    +coeff_weno5_alpha[2][2]*f4;
    
    alpha0=coeff_weno5_c[0]/pow((1.0e-6+is0),2);
    alpha1=coeff_weno5_c[1]/pow((1.0e-6+is1),2);
    alpha2=coeff_weno5_c[2]/pow((1.0e-6+is2),2);
    
    
    alphaSum=alpha0+alpha1+alpha2;
    
    omega0=alpha0/alphaSum;
    omega1=alpha1/alphaSum;
    omega2=alpha2/alphaSum;
    
    *fout=omega0*q30+omega1*q31+omega2*q32;
    
}

