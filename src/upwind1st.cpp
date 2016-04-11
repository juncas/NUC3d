//
//  upwind1st.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/12/12.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#include "upwind1st.hpp"
#include "schemes.hpp"
nuc3d::upwind1st::upwind1st(){
    
}

nuc3d::upwind1st::~upwind1st()
{
    
}


void nuc3d::upwind1st::interpolationInner(const Field & fieldIN,
                                          const int uw,//uw=(1,-1)
                                          Field & fieldOUT,
                                          const int tilesize)
{
    switch(uw)
    {
        case 1:
            upwind1stp(fieldIN, fieldOUT, tilesize);
            break;
        case -1:
            upwind1stn(fieldIN, fieldOUT, tilesize);
            break;
        default:
            std::cout<<"weno error: no such direction"<<std::endl;
            exit(-1);
    }
    
}

void nuc3d::upwind1st::upwind1stp(const Field & fieldIN,
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
            for(int i=ibeg;i<iend;i++)
            {
                int idx_f=nx*ny*k+nx*j+i;
                int idx_rf=nx0*ny0*k+nx0*j+i+1;
                
                pOut[idx_rf]=pIn[idx_f];
                
            }
        }
    }
}

void nuc3d::upwind1st::upwind1stn(const Field & fieldIN,
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
            for(int i=ibeg;i<iend;i++)
            {
                int idx_f=nx*ny*k+nx*j+i;
                int idx_rf=nx0*ny0*k+nx0*j+i+1;
                
                pOut[idx_rf]=pIn[idx_f+1];
                
            }
        }
    }
}

void nuc3d::upwind1st::interpolationBoundaryL(const Field & fieldIN,
                                              const Field & boundaryL,
                                              const int uw,//uw=(1,-1)
                                              Field & fieldOUT,
                                              const int tilesize)
{
    switch(uw)
    {
        case 1:
            upwind1stpBL(fieldIN,boundaryL, fieldOUT, tilesize);
            break;
        case -1:
            upwind1stnBL(fieldIN,boundaryL, fieldOUT, tilesize);
            break;
        default:
            std::cout<<"weno error: no such direction"<<std::endl;
            exit(-1);
    }
    
}

void nuc3d::upwind1st::upwind1stpBL(const Field & fieldIN,
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
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i;
                if((i-1)<0)
                {
                    int idx_BND=nxBND*nyBND*k+nxBND*j+(i-1+nxBND);
                    
                    pOut[idx_rf]=pBND[idx_BND];
                }
                else
                {
                    int idx_f=nx*ny*k+nx*j+i-1;
                    pOut[idx_rf]=pIn[idx_f];
                }
            }
        }
    }
}

void nuc3d::upwind1st::upwind1stnBL(const Field & fieldIN,
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
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i;
                if(i<0)
                {
                    int idx_BND=nxBND*nyBND*k+nxBND*j+i+nxBND;
                    
                    pOut[idx_rf]=pBND[idx_BND];
                }
                else
                {
                    int idx_f=nx*ny*k+nx*j+i;
                    pOut[idx_rf]=pIn[idx_f];
                }
            }
        }
    }
}

void nuc3d::upwind1st::interpolationBoundaryR(const Field & fieldIN,
                                              const Field & boundaryR,
                                              const int uw,//uw=(1,-1)
                                              Field & fieldOUT,
                                              const int tilesize)
{
    switch(uw)
    {
        case 1:
            upwind1stpBR(fieldIN,boundaryR, fieldOUT, tilesize);
            break;
        case -1:
            upwind1stnBR(fieldIN,boundaryR, fieldOUT, tilesize);
            break;
        default:
            std::cout<<"weno error: no such direction"<<std::endl;
            exit(-1);
    }
    
}


void nuc3d::upwind1st::upwind1stpBR(const Field & fieldIN,
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
    
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i+1;
                
                if(i>=nx)
                {
                    int idx_BND=nxBND*nyBND*k+nxBND*j+i-nx;
                    
                    pOut[idx_rf]=pBND[idx_BND];
                }
                else
                {
                    int idx_f=nx*ny*k+nx*j+i;
                    pOut[idx_rf]=pIn[idx_f];
                }
                
            }
        }
    }
}

void nuc3d::upwind1st::upwind1stnBR(const Field & fieldIN,
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
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i+1;
                if((i+1)>=nx)
                {
                    int idx_BND=nxBND*nyBND*k+nxBND*j+i+1-nx;
                    
                    pOut[idx_rf]=pBND[idx_BND];
                }
                else
                {
                    int idx_f=nx*ny*k+nx*j+i+1;
                    pOut[idx_rf]=pIn[idx_f];
                }
                
            }
        }
    }
}
