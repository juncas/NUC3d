//
//  weno3js.cpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/13.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//

#include "weno3js.hpp"

#include "schemes.hpp"

nuc3d::weno3js::weno3js():
ss(1.0e-6),
p(2)
{
    
}

nuc3d::weno3js::~weno3js()
{
    
}

void nuc3d::weno3js::interpolationInner(const Field & fieldIN,
                                        const int uw,//uw=(1,-1)
                                        Field & fieldOUT,
                                        const int tilesize)
{
    switch(uw)
    {
        case 1:
            weno3jsp(fieldIN, fieldOUT, tilesize);
            break;
        case -1:
            weno3jsn(fieldIN, fieldOUT, tilesize);
            break;
        default:
            std::cout<<"weno error: no such direction"<<std::endl;
            exit(-1);
    }
    
}

void nuc3d::weno3js::weno3jsp(const Field & fieldIN,
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
                
                double is0=std::pow((pIn[idx_f-1]-pIn[idx_f]),2);
                double is1=std::pow((pIn[idx_f]-pIn[idx_f+1]),2);
                
                double q20=-0.5*pIn[idx_f-1]+1.5*pIn[idx_f];
                double q21=0.5*pIn[idx_f]+0.5*pIn[idx_f+1];
                
                double aa0=1.0/3.0/std::pow(ss+is0,p);
                double aa1=2.0/3.0/std::pow(ss+is1,p);
                
                double w0=aa0/(aa0+aa1);
                double w1=aa1/(aa0+aa1);
                
                pOut[idx_rf]=w0*q20+w1*q21;
            }
        }
    }
}

void nuc3d::weno3js::weno3jsn(const Field & fieldIN,
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
                
                double is0=std::pow((pIn[idx_f+2]-pIn[idx_f+1]),2);
                double is1=std::pow((pIn[idx_f+1]-pIn[idx_f]),2);
                
                double q20=-0.5*pIn[idx_f+2]+1.5*pIn[idx_f+1];
                double q21=0.5*pIn[idx_f+1]+0.5*pIn[idx_f];
                
                double aa0=1.0/3.0/std::pow(ss+is0,p);
                double aa1=2.0/3.0/std::pow(ss+is1,p);
                
                double w0=aa0/(aa0+aa1);
                double w1=aa1/(aa0+aa1);
                
                pOut[idx_rf]=w0*q20+w1*q21;

            }
        }
    }
}

void nuc3d::weno3js::interpolationBoundaryL(const Field & fieldIN,
                                            const Field & boundaryL,
                                            const int uw,//uw=(1,-1)
                                            Field & fieldOUT,
                                            const int tilesize)
{
    switch(uw)
    {
        case 1:
            weno3jspBL(fieldIN,boundaryL, fieldOUT, tilesize);
            break;
        case -1:
            weno3jsnBL(fieldIN,boundaryL, fieldOUT, tilesize);
            break;
        default:
            std::cout<<"weno error: no such direction"<<std::endl;
            exit(-1);
    }
    
}

void nuc3d::weno3js::weno3jspBL(const Field & fieldIN,
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
    
    double f[3];
    
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i;
                for(int z=-1;z<=1;z++)
                {
                    if((i+z-1)<0)
                    {
                        int idx_BND=nxBND*nyBND*k+nxBND*j+(i+z-1+nxBND);
                        
                        f[1+z]=pBND[idx_BND];
                    }
                    else
                    {
                        int idx_f=nx*ny*k+nx*j+i+z-1;
                        f[1+z]=pIn[idx_f];
                    }
                }
                
                
                double is0=std::pow((f[0]-f[1]),2);
                double is1=std::pow((f[1]-f[2]),2);

                double q20=-0.5*f[0]+1.5*f[1];
                double q21=0.5*f[1]+0.5*f[2];
                
                double aa0=1.0/3.0/std::pow(ss+is0,p);
                double aa1=2.0/3.0/std::pow(ss+is1,p);
                
                double w0=aa0/(aa0+aa1);
                double w1=aa1/(aa0+aa1);
                
                pOut[idx_rf]=w0*q20+w1*q21;
            }
        }
    }
}

void nuc3d::weno3js::weno3jsnBL(const Field & fieldIN,
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
    
    double f[3];
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i;
                
                for(int z=-1;z<=1;z++)
                {
                    if((i-z)<0)
                    {
                        int idx_BND=nxBND*nyBND*k+nxBND*j+(i-z+nxBND);
                        
                        f[1+z]=pBND[idx_BND];
                    }
                    else
                    {
                        int idx_f=nx*ny*k+nx*j+i-z;
                        f[1+z]=pIn[idx_f];
                    }
                }
                
                double is0=std::pow((f[0]-f[1]),2);
                double is1=std::pow((f[1]-f[2]),2);

                 
                double q20=-0.5*f[0]+1.5*f[1];
                double q21=0.5*f[1]+0.5*f[2];
                
                double aa0=1.0/3.0/std::pow(ss+is0,p);
                double aa1=2.0/3.0/std::pow(ss+is1,p);
                
                double w0=aa0/(aa0+aa1);
                double w1=aa1/(aa0+aa1);
                
                pOut[idx_rf]=w0*q20+w1*q21;
            }
        }
    }
}

void nuc3d::weno3js::interpolationBoundaryR(const Field & fieldIN,
                                            const Field & boundaryR,
                                            const int uw,//uw=(1,-1)
                                            Field & fieldOUT,
                                            const int tilesize)
{
    switch(uw)
    {
        case 1:
            weno3jspBR(fieldIN,boundaryR, fieldOUT, tilesize);
            break;
        case -1:
            weno3jsnBR(fieldIN,boundaryR, fieldOUT, tilesize);
            break;
        default:
            std::cout<<"weno error: no such direction"<<std::endl;
            exit(-1);
    }
    
}


void nuc3d::weno3js::weno3jspBR(const Field & fieldIN,
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
    
    double f[3];
    
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i+1;
                for(int z=-1;z<=1;z++)
                {
                    if((i+z)>=nx)
                    {
                        int idx_BND=nxBND*nyBND*k+nxBND*j+i+z-nx;
                        
                        f[1+z]=pBND[idx_BND];
                    }
                    else
                    {
                        int idx_f=nx*ny*k+nx*j+i+z;
                        f[1+z]=pIn[idx_f];
                    }
                }
                
                double is0=std::pow((f[0]-f[1]),2);
                double is1=std::pow((f[1]-f[2]),2);

                double q20=-0.5*f[0]+1.5*f[1];
                double q21=0.5*f[1]+0.5*f[2];
                
                double aa0=1.0/3.0/std::pow(ss+is0,p);
                double aa1=2.0/3.0/std::pow(ss+is1,p);
                
                double w0=aa0/(aa0+aa1);
                double w1=aa1/(aa0+aa1);
                
                pOut[idx_rf]=w0*q20+w1*q21;
            }
        }
    }
}

void nuc3d::weno3js::weno3jsnBR(const Field & fieldIN,
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
    
    double f[3];
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i+1;
                
                for(int z=-1;z<=1;z++)
                {
                    if((i-z+1)>=nx)
                    {
                        int idx_BND=nxBND*nyBND*k+nxBND*j+i-z+1-nx;
                        
                        f[1+z]=pBND[idx_BND];
                    }
                    else
                    {
                        int idx_f=nx*ny*k+nx*j+i-z+1;
                        f[1+z]=pIn[idx_f];
                    }
                }
                
                double is0=std::pow((f[0]-f[1]),2);
                double is1=std::pow((f[1]-f[2]),2);

                double q20=-0.5*f[0]+1.5*f[1];
                double q21=0.5*f[1]+0.5*f[2];
                
                double aa0=1.0/3.0/std::pow(ss+is0,p);
                double aa1=2.0/3.0/std::pow(ss+is1,p);
                
                double w0=aa0/(aa0+aa1);
                double w1=aa1/(aa0+aa1);
                
                pOut[idx_rf]=w0*q20+w1*q21;
            }
        }
    }
}

