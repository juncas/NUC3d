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
                                          const int dim0,
                                          const int dim1,
                                          const int dim2,
                                          const int uw,//uw=(1,-1)
                                          Field & fieldOUT,
                                          const int tilesize)
{
    double h;
    
    
    int nx;
    int ny;
    int nz;
    
    nx=fieldIN.getSizeX();
    ny=fieldIN.getSizeY();
    nz=fieldIN.getSizeZ();
    
    int ibeg=(tilesize-1)*dim0;
    int iend=nx-tilesize*dim0;
    int jbeg=(tilesize-1)*dim1;
    int jend=ny-tilesize*dim1;
    int kbeg=(tilesize-1)*dim2;
    int kend=nz-tilesize*dim2;
    
    if( ((nx+dim0)==(fieldOUT.getSizeX()))&&
       ((ny+dim1)==(fieldOUT.getSizeY()))&&
       ((nz+dim2)==(fieldOUT.getSizeZ())) )
    {
        for(int k=kbeg;k<kend;k++)
        {
            for(int j=jbeg;j<jend;j++)
            {
                for(int i=ibeg;i<iend;i++)
                {
                    double hl=fieldIN.getValue(i,j,k);
                    double hr=fieldIN.getValue(i+dim0,j+dim1,k+dim2);
                    
                    h=0.5*(hl+hr);
                    
                    fieldOUT.setValue(i+dim0,j+dim1,k+dim2,h);
                }
            }
        }
    }
}


void nuc3d::upwind1st::interpolationBoundaryL(const Field & fieldIN,
                                              const Field & boundaryL,
                                              const int dim0,
                                              const int dim1,
                                              const int dim2,
                                              const int uw,//uw=(1,-1)
                                              Field & fieldOUT,
                                              const int tilesize)
{
    double h;
    double hl;
    double hr;
    
    int nx;
    int ny;
    int nz;
    
    nx=fieldIN.getSizeX();
    ny=fieldIN.getSizeY();
    nz=fieldIN.getSizeZ();
    
    int ibeg=0;
    int iend=dim0*(tilesize-nx)+nx;
    int jbeg=0;
    int jend=dim1*(tilesize-ny)+ny;
    int kbeg=0;
    int kend=dim2*(tilesize-nz)+nz;
    
    
    if( ((nx+dim0)==(fieldOUT.getSizeX()))&&
       ((ny+dim1)==(fieldOUT.getSizeY()))&&
       ((nz+dim2)==(fieldOUT.getSizeZ())) )
    {
        for(int k=kbeg;k<kend;k++)
        {
            for(int j=jbeg;j<jend;j++)
            {
                for(int i=ibeg;i<iend;i++)
                {
                    int itemp0=i-dim0;
                    int jtemp0=j-dim1;
                    int ktemp0=k-dim2;
                    
                    if(itemp0<0||jtemp0<0||ktemp0<0)
                    {
                        int itemp1=itemp0+dim0*boundaryL.getSizeX();
                        int jtemp1=jtemp0+dim1*boundaryL.getSizeY();
                        int ktemp1=ktemp0+dim2*boundaryL.getSizeZ();
                        
                        hl=boundaryL.getValue(itemp1,jtemp1,ktemp1);
                    }
                    else
                    {
                        hl=fieldIN.getValue(itemp0,jtemp0,ktemp0);
                    }
                    
                    hr=fieldIN.getValue(i, j, k);
                    
                    h=(hl+hr)*0.5;
                    
                    fieldOUT.setValue(i,j,k,h);
                }
            }
        }
    }
}

void nuc3d::upwind1st::interpolationBoundaryR(const Field & fieldIN,
                                              const Field & boundaryR,
                                              const int dim0,
                                              const int dim1,
                                              const int dim2,
                                              const int uw,//uw=(1,-1)
                                              Field & fieldOUT,
                                              const int tilesize)
{
    double h;
    double hl,hr;
    
    
    int nx;
    int ny;
    int nz;
    
    nx=fieldIN.getSizeX();
    ny=fieldIN.getSizeY();
    nz=fieldIN.getSizeZ();
    
    int ibeg=dim0*(nx-tilesize);
    int iend=nx;
    int jbeg=dim1*(ny-tilesize);
    int jend=ny;
    int kbeg=dim2*(nz-tilesize);
    int kend=nz;
    
    if( ((nx+dim0)==(fieldOUT.getSizeX()))&&
       ((ny+dim1)==(fieldOUT.getSizeY()))&&
       ((nz+dim2)==(fieldOUT.getSizeZ())) )
    {
        for(int k=kbeg;k<kend;k++)
        {
            for(int j=jbeg;j<jend;j++)
            {
                for(int i=ibeg;i<iend;i++)
                {
                    int itemp0=i+dim0;
                    int jtemp0=j+dim1;
                    int ktemp0=k+dim2;
                    
                    if(itemp0>=nx||jtemp0>=ny||ktemp0>=nz)
                    {
                        int itemp1=itemp0-dim0*nx;
                        int jtemp1=jtemp0-dim1*ny;
                        int ktemp1=ktemp0-dim2*nz;
                        
                        hr=boundaryR.getValue(itemp1,jtemp1,ktemp1);
                    }
                    else
                    {
                        hr=fieldIN.getValue(itemp0,jtemp0,ktemp0);
                    }
                    
                    hl=fieldIN.getValue(i,j,k);
                    
                    h=(hl+hr)*0.5;
                    
                    fieldOUT.setValue(i+dim0,j+dim1,k+dim2,h);
                }
            }
        }
    }
}
