#ifndef weno5js_cpp
#define weno5js_cpp
#include "weno5js.h"
nuc3d::weno5js::weno5js():
ss(1.0e-6),
p(2)
{

}

nuc3d::weno5js::~weno5js()
{

}

double nuc3d::weno5js::weno5jsInterpolation(const double *f)
{
    double omega0,omega1,omega2;
    double alpha0,alpha1,alpha2,alphaSum;
    double is0,is1,is2;
    double q30,q31,q32;
    double tau;

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
    
    alpha0=coeff_weno5_c[0]/pow((ss+is0),p);
    alpha1=coeff_weno5_c[1]/pow((ss+is1),p);
    alpha2=coeff_weno5_c[2]/pow((ss+is2),p);
    
    tau=std::abs(is2-is0);
    
    alpha0=coeff_weno5_c[0]*(1.0+std::pow(tau/(is0+ss),p));
    alpha1=coeff_weno5_c[1]*(1.0+std::pow(tau/(is1+ss),p));
    alpha2=coeff_weno5_c[2]*(1.0+std::pow(tau/(is2+ss),p));
    
    alphaSum=alpha0+alpha1+alpha2;

    omega0=alpha0/alphaSum;
    omega1=alpha1/alphaSum;
    omega2=alpha2/alphaSum;

    return omega0*q30+omega1*q31+omega2*q32;
    
}

void nuc3d::weno5js::interpolationInner(
  const Field & fieldIN,
  const int dim0,
  const int dim1,
  const int dim2,
  const int uw,//uw=(1,-1)
        Field & fieldOUT
)
{
    double flux[5];
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
                        for(int z=-2;z<=2;z++)
                        {
                            int stride=(1-uw)/2+z*uw;
                            int itemp=i+stride*dim0;
                            int jtemp=j+stride*dim1;
                            int ktemp=k+stride*dim2;

                            flux[2+z]=fieldIN.getValue(itemp,jtemp,ktemp);
                        }
                        h=weno5jsInterpolation(flux);
                        
                        fieldOUT.setValue(i+dim0,j+dim1,k+dim2,h);
                    }
                    
                }
            }
        }
    }
}


void nuc3d::weno5js::interpolationBoundary(
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
    double flux[5];
    double h;


    int nx;
    int ny;
    int nz;

    nx=fieldIN.getSizeX();
    ny=fieldIN.getSizeY();
    nz=fieldIN.getSizeZ();
    
    int strideL=3-(1-uw)/2;
    int strideR=1+(1-uw)/2;

    
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
                    
                        h=weno5jsInterpolation(flux);
    
                        fieldOUT.setValue(i,j,k,h);

                    }
                    else if(!(((i+3*dim0)<nx)&&((j+3*dim1)<ny)&&((k+3*dim2)<nz)))
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
                    
                    
                        h=weno5jsInterpolation(flux);
    
                        fieldOUT.setValue(i+dim0,j+dim1,k+dim2,h);
                    }
                }
            }
        }
    }
}

#endif
