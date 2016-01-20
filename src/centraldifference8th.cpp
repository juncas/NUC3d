//
//  centraldifference8th.cpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/20.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//

#include "centraldifference8th.hpp"
#include "schemes.hpp"

nuc3d::centraldifference8th::centraldifference8th()
{
    
}

nuc3d::centraldifference8th::~centraldifference8th()
{
    
}

void nuc3d::centraldifference8th::differentialInner(const Field & fieldIN,
                                                    Field & fieldOUT,
                                                    const int tilesize)
{
    
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    
    int nx=fieldIN.getSizeX();
    int ny=fieldIN.getSizeY();
    int nz=fieldIN.getSizeZ();
    
    int ibeg=tilesize;
    int iend=nx-tilesize;
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
                int idx_f=nx*ny*k+nx*j+i;
                
                pOut[idx_f]=coeff_cd8_alpha[3]*(pIn[idx_f+4]-pIn[idx_f-4])
                +coeff_cd8_alpha[2]*(pIn[idx_f+3]-pIn[idx_f-3])
                +coeff_cd8_alpha[1]*(pIn[idx_f+2]-pIn[idx_f-2])
                +coeff_cd8_alpha[0]*(pIn[idx_f+1]-pIn[idx_f-1]);
            }
        }
    }
}


void nuc3d::centraldifference8th::differentialBoundaryL(const Field & fieldIN,
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
    
    int nxBND=boundaryL.getSizeX();
    int nyBND=boundaryL.getSizeY();
    int nzBND=boundaryL.getSizeZ();
    
    int ibeg=0;
    int iend=tilesize;
    int jbeg=0;
    int jend=ny;
    int kbeg=0;
    int kend=nz;
    
    double f[8];
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                for(int z=-4;z<0;z++)
                {
                    if((i+z)<0)
                    {
                        int idx_BND=nxBND*nyBND*k+nxBND*j+(i+z+nxBND);
                        
                        f[4+z]=pBND[idx_BND];
                    }
                    else
                    {
                        int idx_f=nx*ny*k+nx*j+i+z;
                        
                        f[4+z]=pIn[idx_f];
                    }
                }
                
                for(int z=1;z<=4;z++)
                {
                    int idx_f=nx*ny*k+nx*j+i+z;
                    
                    f[3+z]=pIn[idx_f];
                }
                
                int idx_f_out=nx*ny*k+nx*j+i;
                
                pOut[idx_f_out]=coeff_cd8_alpha[3]*(f[7]-f[0])
                +coeff_cd8_alpha[2]*(f[6]-f[1])
                +coeff_cd8_alpha[1]*(f[5]-f[2])
                +coeff_cd8_alpha[0]*(f[4]-f[3]);
                
            }
        }
    }
}


void nuc3d::centraldifference8th::differentialBoundaryR(const Field & fieldIN,
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
    
    int nxBND=boundaryR.getSizeX();
    int nyBND=boundaryR.getSizeY();
    int nzBND=boundaryR.getSizeZ();
    
    int ibeg=nx-tilesize;
    int iend=nx;
    int jbeg=0;
    int jend=ny;
    int kbeg=0;
    int kend=nz;
    
    double f[8];
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                for(int z=-4;z<0;z++)
                {
                    int idx_f=nx*ny*k+nx*j+i+z;
                    
                    f[4+z]=pIn[idx_f];
                }
                
                for(int z=1;z<=4;z++)
                {
                    if((i+z)>=nx)
                    {
                        int idx_BND=nxBND*nyBND*k+nxBND*j+(i+z-nx);
                        
                        f[3+z]=pBND[idx_BND];
                    }
                    else
                    {
                        int idx_f=nx*ny*k+nx*j+i+z;
                        
                        f[3+z]=pIn[idx_f];
                    }
                }
                
                int idx_f_out=nx*ny*k+nx*j+i;
                
                pOut[idx_f_out]=coeff_cd8_alpha[3]*(f[7]-f[0])
                +coeff_cd8_alpha[2]*(f[6]-f[1])
                +coeff_cd8_alpha[1]*(f[5]-f[2])
                +coeff_cd8_alpha[0]*(f[4]-f[3]);
                
            }
        }
    }
}
