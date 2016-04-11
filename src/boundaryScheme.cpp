//
//  boundaryScheme.cpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/9.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//

#include "boundaryScheme.hpp"

#include "schemes.hpp"

nuc3d::boundaryScheme::boundaryScheme()
{
    
}

nuc3d::boundaryScheme::~boundaryScheme()
{
    
}

void nuc3d::boundaryScheme::differentialInner(const Field & fieldIN,
                                              Field & fieldOUT,
                                              const int tilesize)
{
    
}


void nuc3d::boundaryScheme::differentialBoundaryL(const Field & fieldIN,
                                                  const Field & boundaryL,
                                                  Field & fieldOUT,
                                                  const int tilesize)
{
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    
    int nx=fieldIN.getSizeX();
    int ny=fieldIN.getSizeY();
    int nz=fieldIN.getSizeZ();
    
    int ibeg=0;
    int iend=tilesize;
    int jbeg=0;
    int jend=ny;
    int kbeg=0;
    int kend=nz;
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_f_out=nx*ny*k+nx*j+i;
                
                (this->*myBCschemeL[i])(pIn,pOut,idx_f_out);
                
            }
        }
    }
}


void nuc3d::boundaryScheme::differentialBoundaryR(const Field & fieldIN,
                                                  const Field & boundaryR,
                                                  Field & fieldOUT,
                                                  const int tilesize)
{
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    
    int nx=fieldIN.getSizeX();
    int ny=fieldIN.getSizeY();
    int nz=fieldIN.getSizeZ();
    
    int ibeg=nx-tilesize;
    int iend=nx;
    int jbeg=0;
    int jend=ny;
    int kbeg=0;
    int kend=nz;
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_f_out=nx*ny*k+nx*j+i;
                
                (this->*myBCschemeR[iend-i-1])(pIn,pOut,idx_f_out);
                
            }
        }
    }
}

void nuc3d::boundaryScheme::firstOrderL(double *f_in,double *f_out,int &idx0)
{
    f_out[idx0]=f_in[idx0+1]-f_in[idx0];
}

void nuc3d::boundaryScheme::firstOrderR(double *f_in,double *f_out,int &idx0)
{
    f_out[idx0]=f_in[idx0]-f_in[idx0-1];
}

void nuc3d::boundaryScheme::secondOrder(double *f_in,double *f_out,int &idx0)
{
    f_out[idx0]=coeff_cd2_alpha*(f_in[idx0+1]-f_in[idx0-1]);
}

void nuc3d::boundaryScheme::fourthOrder(double *f_in,double *f_out,int &idx0)
{
    f_out[idx0]=coeff_cd4_alpha[1]*(f_in[idx0+2]-f_in[idx0-2])
    +coeff_cd4_alpha[0]*(f_in[idx0+1]-f_in[idx0-1]);
}

void nuc3d::boundaryScheme::sixthOrder(double *f_in,double *f_out,int &idx0)
{
    f_out[idx0]=coeff_cd6_alpha[2]*(f_in[idx0+3]-f_in[idx0-3])
    +coeff_cd6_alpha[1]*(f_in[idx0+2]-f_in[idx0-2])
    +coeff_cd4_alpha[0]*(f_in[idx0+1]-f_in[idx0-1]);

}
