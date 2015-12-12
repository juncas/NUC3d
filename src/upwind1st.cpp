//
//  upwind1st.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/12/12.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#include "upwind1st.hpp"
nuc3d::upwind1st::upwind1st(){
    
}

nuc3d::upwind1st::~upwind1st()
{
    
}


void nuc3d::upwind1st::interpolationInner(                                        const Field & fieldIN,
                                        const int dim0,
                                        const int dim1,
                                        const int dim2,
                                        const int uw,//uw=(1,-1)
                                        Field & fieldOUT
                                        )
{
    double h;
    
    
    int nx;
    int ny;
    int nz;
    
    nx=fieldIN.getSizeX();
    ny=fieldIN.getSizeY();
    nz=fieldIN.getSizeZ();
    
    int strideL=2-(1-uw)/2;
    int strideR=2+(1-uw)/2;
    
    if( ((nx+dim0)==(fieldOUT.getSizeX()))&&
       ((ny+dim1)==(fieldOUT.getSizeY()))&&
       ((nz+dim2)==(fieldOUT.getSizeZ())) )
    {
        for(int k=0;k<nz;k++)
        {
            for(int j=0;j<ny;j++)
            {
                for(int i=0;i<nx;i++)
                {
                    if((((i-strideL*dim0)>=0)&&((j-strideL*dim1)>=0)&&((k-strideL*dim2)>=0))&&
                       (((i+strideR*dim0)<nx)&&((j+strideR*dim1)<ny)&&((k+strideR*dim2)<nz)))
                    {
                            int stride=(1-uw)/2;
                            int itemp=i+stride*dim0;
                            int jtemp=j+stride*dim1;
                            int ktemp=k+stride*dim2;
                        
                        h=fieldIN.getValue(itemp,jtemp,ktemp);
                        fieldOUT.setValue(i+dim0,j+dim1,k+dim2,h);
                    }
                    
                }
            }
        }
    }
}


void nuc3d::upwind1st::interpolationBoundary(
                                           const Field & fieldIN,
                                           const Field & boundaryL,
                                           const Field & boundaryR,
                                           const int dim0,
                                           const int dim1,
                                           const int dim2,
                                           const int uw,//uw=(1,-1)
                                           Field & fieldOUT
                                           )
{
    double h;
    
    
    int nx;
    int ny;
    int nz;
    
    nx=fieldIN.getSizeX();
    ny=fieldIN.getSizeY();
    nz=fieldIN.getSizeZ();
        
    if( ((nx+dim0)==(fieldOUT.getSizeX()))&&
       ((ny+dim1)==(fieldOUT.getSizeY()))&&
       ((nz+dim2)==(fieldOUT.getSizeZ())) )
    {
        for(int k=0;k<nz;++k)
        {
            for(int j=0;j<ny;++j)
            {
                for(int i=0;i<nx;++i)
                {
                    if(!(((i-3*dim0)>=0)&&((j-3*dim1)>=0)&&((k-3*dim2)>=0)))
                    {
                        
                            int stride=(1-uw)/2-1;
                            int itemp0=i+stride*dim0;
                            int jtemp0=j+stride*dim1;
                            int ktemp0=k+stride*dim2;
                            
                            if(itemp0<0||jtemp0<0||ktemp0<0)
                            {
                                int itemp1=itemp0+dim0*boundaryL.getSizeX();
                                int jtemp1=jtemp0+dim1*boundaryL.getSizeY();
                                int ktemp1=ktemp0+dim2*boundaryL.getSizeZ();
                                
                                h=boundaryL.getValue(itemp1,jtemp1,ktemp1);
                            }
                            else
                            {
                                h=fieldIN.getValue(itemp0,jtemp0,ktemp0);
                            }
                        
                        fieldOUT.setValue(i,j,k,h);
                        
                    }
                    else if(!(((i+3*dim0)<nx)&&((j+3*dim1)<ny)&&((k+3*dim2)<nz)))
                    {
                        
                            int stride=(1-uw)/2;
                            int itemp0=i+stride*dim0;
                            int jtemp0=j+stride*dim1;
                            int ktemp0=k+stride*dim2;
                            
                            if(itemp0>=nx||jtemp0>=ny||ktemp0>=nz)
                            {
                                int itemp1=itemp0-dim0*nx;
                                int jtemp1=jtemp0-dim1*ny;
                                int ktemp1=ktemp0-dim2*nz;
                                
                                h=boundaryR.getValue(itemp1,jtemp1,ktemp1);
                            }
                            else
                            {
                                h=fieldIN.getValue(itemp0,jtemp0,ktemp0);
                            }
                        
                        fieldOUT.setValue(i+dim0,j+dim1,k+dim2,h);
                    }
                }
            }
        }
    }
}
