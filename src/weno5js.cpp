#ifndef weno5js_cpp
#define weno5js_cpp
#include "weno5js.h"
#include "schemes.hpp"

nuc3d::weno5js::weno5js():
ss(1.0e-6),
p(2)
{
    
}

nuc3d::weno5js::~weno5js()
{
    
}

void nuc3d::weno5js::interpolationInner(const Field & fieldIN,
                                        const int dim0,
                                        const int dim1,
                                        const int dim2,
                                        const int uw,//uw=(1,-1)
                                        Field & fieldOUT,
                                        const int tilesize)
{
    double flux[5];
    double h;
    const double *pIn = nullptr;
    pIn=fieldIN.getDataPtr();
    
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
                    for(int z=-2;z<=2;z++)
                    {
                        int stride=(1-uw)/2+z*uw;
                        int itemp=i+stride*dim0;
                        int jtemp=j+stride*dim1;
                        int ktemp=k+stride*dim2;
                        int idx=ktemp*nx*ny+jtemp*nx+itemp;
                        
                        //flux[2+z]=fieldIN.getValue(itemp,jtemp,ktemp);
                        flux[2+z]=pIn[idx];
                    }
                    h=nuc3d::weno5jsInterpolation(flux,ss,p);
                    
                    fieldOUT.setValue(i+dim0,j+dim1,k+dim2,h);
                }
            }
        }
    }
    else
    {
        std::cout<<"Error size in interpolationInner"<<dim0<<" "<<dim1<<" "<<dim2<<std::endl;
    }
}


void nuc3d::weno5js::interpolationBoundaryL(const Field & fieldIN,
                                            const Field & boundaryL,
                                            const int dim0,
                                            const int dim1,
                                            const int dim2,
                                            const int uw,//uw=(1,-1)
                                            Field & fieldOUT,
                                            const int tilesize)
{
    double flux[5];
    double h;
    
    
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
                    for(int z=-2;z<=2;z++)
                    {
                        int stride=(1-uw)/2+z*uw-1;
                        int itemp0=i+stride*dim0;
                        int jtemp0=j+stride*dim1;
                        int ktemp0=k+stride*dim2;
                        
                        if(itemp0<0||jtemp0<0||ktemp0<0)
                        {
                            int itemp1=itemp0+dim0*boundaryL.getSizeX();
                            int jtemp1=jtemp0+dim1*boundaryL.getSizeY();
                            int ktemp1=ktemp0+dim2*boundaryL.getSizeZ();
                            
                            flux[2+z]=boundaryL.getValue(itemp1,jtemp1,ktemp1);
                        }
                        else
                        {
                            flux[2+z]=fieldIN.getValue(itemp0,jtemp0,ktemp0);
                        }
                        
                    }
                    
                    h=nuc3d::weno5jsInterpolation(flux,ss,p);
                    
                    fieldOUT.setValue(i,j,k,h);
                }
            }
        }
    }
}

void nuc3d::weno5js::interpolationBoundaryR(const Field & fieldIN,
                                            const Field & boundaryR,
                                            const int dim0,
                                            const int dim1,
                                            const int dim2,
                                            const int uw,//uw=(1,-1)
                                            Field & fieldOUT,
                                            const int tilesize)
{
    double flux[5];
    double h;
    
    
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
                    for(int z=-2;z<=2;z++)
                    {
                        int stride=(1-uw)/2+z*uw;
                        int itemp0=i+stride*dim0;
                        int jtemp0=j+stride*dim1;
                        int ktemp0=k+stride*dim2;
                        
                        if(itemp0>=nx||jtemp0>=ny||ktemp0>=nz)
                        {
                            int itemp1=itemp0-dim0*nx;
                            int jtemp1=jtemp0-dim1*ny;
                            int ktemp1=ktemp0-dim2*nz;
                            
                            flux[2+z]=boundaryR.getValue(itemp1,jtemp1,ktemp1);
                        }
                        else
                        {
                            flux[2+z]=fieldIN.getValue(itemp0,jtemp0,ktemp0);
                        }
                        
                    }
                    
                    
                    h=nuc3d::weno5jsInterpolation(flux,ss,p);
                    
                    fieldOUT.setValue(i+dim0,j+dim1,k+dim2,h);
                }
            }
        }
    }
}

#endif
