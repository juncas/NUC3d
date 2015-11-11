#ifndef centraldifference6th_cpp
#define centraldifference6th_cpp
#include "centraldifference6th.h"

nuc3d::centraldifference6th::centraldifference6th()
{

}

nuc3d::centraldifference6th::~centraldifference6th()
{

}

double nuc3d::centraldifference6th::cd6differential(const double *f)
{
    double df;

    df = coeff_cd6_alpha[2]*(f[5]-f[0])
        +coeff_cd6_alpha[1]*(f[4]-f[1])
        +coeff_cd6_alpha[0]*(f[3]-f[2]);
        
    return df;
}

void nuc3d::centraldifference6th::differentialInner(
  const Field & fieldIN,
  const int dim0,
  const int dim1,
  const int dim2,
        Field & fieldOUT
)
{
    double flux[6];
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
                    if((((i-3*dim0)>=0)&&((j-3*dim1)>=0)&&((k-3*dim2)>=0))&&
                        (((i+3*dim0)<nx)&&((j+3*dim1)<ny)&&((k+3*dim2)<nz)))
                    {   
                        for(int z=-3;z<=3;z++)
                        {
                            int stride=z;
                            int itemp=i+stride*dim0;
                            int jtemp=j+stride*dim1;
                            int ktemp=k+stride*dim2;

                            flux[2+z]=fieldIN.getValue(itemp,jtemp,ktemp);
                        }
                    }
                    
                    h= cd6differential(flux);
    
                    fieldOUT.setValue(i+dim0,j+dim1,k+dim2,h);
                }
            }
        }
    }
}


void nuc3d::centraldifference6th::differentialBoundary(
  const Field & fieldIN,
  const Field & boundaryL,
  const Field & boundaryR,
  const int dim0,
  const int dim1,
  const int dim2,
        Field & fieldOUT
)
{
    double flux[6];
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
                        for(int z=-3;z<=3;z++)
                        {
                            int stride=z-1;
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
                    
                        h= cd6differential(flux);
    
                        fieldOUT.setValue(i,j,k,h);

                    }
                    else if(!(((i+3*dim0)<nx)&&((j+3*dim1)<ny)&&((k+3*dim2)<nz)))
                    {
                        for(int z=-3;z<=3;z++)
                        {
                            int stride=z;
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
                    
                    
                        h= cd6differential(flux);
    
                        fieldOUT.setValue(i+dim0,j+dim1,k+dim2,h);
                    }
                }
            }
        }
    }
}

#endif
