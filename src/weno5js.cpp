#ifndef weno5js_cpp
#define weno5js_cpp
#include "weno5js.h"
#include "schemes.hpp"

nuc3d::weno5js::weno5js():
ss(1.0e-6),
p(2.0)
{
    
}

nuc3d::weno5js::~weno5js()
{
    
}

void nuc3d::weno5js::interpolationInner(const Field & fieldIN,
                                        const int uw,//uw=(1,-1)
                                        Field & fieldOUT,
                                        const int tilesize)
{
    switch(uw)
    {
        case 1:
            weno5jsp(fieldIN, fieldOUT, tilesize);
            break;
        case -1:
            weno5jsn(fieldIN, fieldOUT, tilesize);
            break;
        default:
            std::cout<<"weno error: no such direction"<<std::endl;
            exit(-1);
    }
    
}

void nuc3d::weno5js::weno5jsp(const Field & fieldIN,
                              Field & fieldOUT,
                              const int tilesize)
{
    
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    
    double omega0,omega1,omega2;
    double alpha0,alpha1,alpha2,alphaSum;
    double is0,is1,is2;
    double q30,q31,q32;
    
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
                
                is0= coeff_weno5_gamma0*std::pow((pIn[idx_f-2]-2.0*pIn[idx_f-1]+pIn[idx_f]),2)
                +coeff_weno5_gamma1*std::pow((pIn[idx_f-2]-4.0*pIn[idx_f-1]+3.0*pIn[idx_f]),2);
                
                is1= coeff_weno5_gamma0*std::pow((pIn[idx_f-1]-2.0*pIn[idx_f]+pIn[idx_f+1]),2)
                +coeff_weno5_gamma1*std::pow((pIn[idx_f-1]-pIn[idx_f+1]),2);
                
                is2= coeff_weno5_gamma0*std::pow((pIn[idx_f]-2.0*pIn[idx_f+1]+pIn[idx_f+2]),2)
                +coeff_weno5_gamma1*std::pow((3.0*pIn[idx_f]-4.0*pIn[idx_f+1]+pIn[idx_f+2]),2);
                
                
                q30= coeff_weno5_alpha[0][0]*pIn[idx_f-2]
                +coeff_weno5_alpha[0][1]*pIn[idx_f-1]
                +coeff_weno5_alpha[0][2]*pIn[idx_f];
                
                q31= coeff_weno5_alpha[1][0]*pIn[idx_f-1]
                +coeff_weno5_alpha[1][1]*pIn[idx_f]
                +coeff_weno5_alpha[1][2]*pIn[idx_f+1];
                
                q32= coeff_weno5_alpha[2][0]*pIn[idx_f]
                +coeff_weno5_alpha[2][1]*pIn[idx_f+1]
                +coeff_weno5_alpha[2][2]*pIn[idx_f+2];
                
                alpha0=coeff_weno5_c[0]/std::pow((ss+is0),p);
                alpha1=coeff_weno5_c[1]/std::pow((ss+is1),p);
                alpha2=coeff_weno5_c[2]/std::pow((ss+is2),p);
                
                
                alphaSum=alpha0+alpha1+alpha2;
                
                omega0=alpha0/alphaSum;
                omega1=alpha1/alphaSum;
                omega2=alpha2/alphaSum;
                
                pOut[idx_rf]=omega0*q30+omega1*q31+omega2*q32;
            }
        }
    }
}

void nuc3d::weno5js::weno5jsn(const Field & fieldIN,
                              Field & fieldOUT,
                              const int tilesize)
{
    
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    
    double omega0,omega1,omega2;
    double alpha0,alpha1,alpha2,alphaSum;
    double is0,is1,is2;
    double q30,q31,q32;
    
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
                
                is0= coeff_weno5_gamma0*std::pow((pIn[idx_f+3]-2.0*pIn[idx_f+2]+pIn[idx_f+1]),2)
                +coeff_weno5_gamma1*std::pow((pIn[idx_f+3]-4.0*pIn[idx_f+2]+3.0*pIn[idx_f+1]),2);
                
                is1= coeff_weno5_gamma0*std::pow((pIn[idx_f+2]-2.0*pIn[idx_f+1]+pIn[idx_f]),2)
                +coeff_weno5_gamma1*std::pow((pIn[idx_f+2]-pIn[idx_f]),2);
                
                is2= coeff_weno5_gamma0*std::pow((pIn[idx_f+1]-2.0*pIn[idx_f]+pIn[idx_f-1]),2)
                +coeff_weno5_gamma1*std::pow((3.0*pIn[idx_f+1]-4.0*pIn[idx_f]+pIn[idx_f-1]),2);
                
                
                q30= coeff_weno5_alpha[0][0]*pIn[idx_f+3]
                +coeff_weno5_alpha[0][1]*pIn[idx_f+2]
                +coeff_weno5_alpha[0][2]*pIn[idx_f+1];
                
                q31= coeff_weno5_alpha[1][0]*pIn[idx_f+2]
                +coeff_weno5_alpha[1][1]*pIn[idx_f+1]
                +coeff_weno5_alpha[1][2]*pIn[idx_f];
                
                q32= coeff_weno5_alpha[2][0]*pIn[idx_f+1]
                +coeff_weno5_alpha[2][1]*pIn[idx_f]
                +coeff_weno5_alpha[2][2]*pIn[idx_f-1];
                
                alpha0=coeff_weno5_c[0]/std::pow((ss+is0),p);
                alpha1=coeff_weno5_c[1]/std::pow((ss+is1),p);
                alpha2=coeff_weno5_c[2]/std::pow((ss+is2),p);
                
                
                alphaSum=alpha0+alpha1+alpha2;
                
                omega0=alpha0/alphaSum;
                omega1=alpha1/alphaSum;
                omega2=alpha2/alphaSum;
                
                pOut[idx_rf]=omega0*q30+omega1*q31+omega2*q32;
            }
        }
    }
}

void nuc3d::weno5js::interpolationBoundaryL(const Field & fieldIN,
                                            const Field & boundaryL,
                                            const int uw,//uw=(1,-1)
                                            Field & fieldOUT,
                                            const int tilesize)
{
    switch(uw)
    {
        case 1:
            weno5jspBL(fieldIN,boundaryL, fieldOUT, tilesize);
            break;
        case -1:
            weno5jsnBL(fieldIN,boundaryL, fieldOUT, tilesize);
            break;
        default:
            std::cout<<"weno error: no such direction"<<std::endl;
            exit(-1);
    }
    
}

void nuc3d::weno5js::weno5jspBL(const Field & fieldIN,
                                const Field & boundaryL,
                                Field & fieldOUT,
                                const int tilesize)
{
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    double *pBND=boundaryL.getDataPtr();
    
    double omega0,omega1,omega2;
    double alpha0,alpha1,alpha2,alphaSum;
    double is0,is1,is2;
    double q30,q31,q32;
    
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
    
    double f[5];
    
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i;
                for(int z=-2;z<=2;z++)
                {
                    if((i+z-1)<0)
                    {
                        int idx_BND=nxBND*nyBND*k+nxBND*j+(i+z-1+nxBND);
                        
                        f[2+z]=pBND[idx_BND];
                    }
                    else
                    {
                        int idx_f=nx*ny*k+nx*j+i+z-1;
                        f[2+z]=pIn[idx_f];
                    }
                }
                
                
                is0= coeff_weno5_gamma0*std::pow((    f[0]-2.0*f[1]+    f[2]),2)
                +coeff_weno5_gamma1*std::pow((    f[0]-4.0*f[1]+3.0*f[2]),2);
                
                is1= coeff_weno5_gamma0*std::pow((    f[1]-2.0*f[2]+    f[3]),2)
                +coeff_weno5_gamma1*std::pow((    f[1]-             f[3]),2);
                
                is2= coeff_weno5_gamma0*std::pow((    f[2]-2.0*f[3]+    f[4]),2)
                +coeff_weno5_gamma1*std::pow((3.0*f[2]-4.0*f[3]+    f[4]),2);
                
                
                q30= coeff_weno5_alpha[0][0]*f[0]
                +coeff_weno5_alpha[0][1]*f[1]
                +coeff_weno5_alpha[0][2]*f[2];
                
                q31= coeff_weno5_alpha[1][0]*f[1]
                +coeff_weno5_alpha[1][1]*f[2]
                +coeff_weno5_alpha[1][2]*f[3];
                
                q32= coeff_weno5_alpha[2][0]*f[2]
                +coeff_weno5_alpha[2][1]*f[3]
                +coeff_weno5_alpha[2][2]*f[4];
                
                alpha0=coeff_weno5_c[0]/std::pow((ss+is0),p);
                alpha1=coeff_weno5_c[1]/std::pow((ss+is1),p);
                alpha2=coeff_weno5_c[2]/std::pow((ss+is2),p);
                
                
                alphaSum=alpha0+alpha1+alpha2;
                
                omega0=alpha0/alphaSum;
                omega1=alpha1/alphaSum;
                omega2=alpha2/alphaSum;
                
                pOut[idx_rf]=omega0*q30+omega1*q31+omega2*q32;
            }
        }
    }
}

void nuc3d::weno5js::weno5jsnBL(const Field & fieldIN,
                                const Field & boundaryL,
                                Field & fieldOUT,
                                const int tilesize)
{
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    double *pBND=boundaryL.getDataPtr();
    
    double omega0,omega1,omega2;
    double alpha0,alpha1,alpha2,alphaSum;
    double is0,is1,is2;
    double q30,q31,q32;
    
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
    
    double f[5];
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i;
                
                for(int z=-2;z<=2;z++)
                {
                    if((i-z)<0)
                    {
                        int idx_BND=nxBND*nyBND*k+nxBND*j+(i-z+nxBND);
                        
                        f[2+z]=pBND[idx_BND];
                    }
                    else
                    {
                        int idx_f=nx*ny*k+nx*j+i-z;
                        f[2+z]=pIn[idx_f];
                    }
                }
                
                
                is0= coeff_weno5_gamma0*std::pow((    f[0]-2.0*f[1]+    f[2]),2)
                +coeff_weno5_gamma1*std::pow((    f[0]-4.0*f[1]+3.0*f[2]),2);
                
                is1= coeff_weno5_gamma0*std::pow((    f[1]-2.0*f[2]+    f[3]),2)
                +coeff_weno5_gamma1*std::pow((    f[1]-             f[3]),2);
                
                is2= coeff_weno5_gamma0*std::pow((    f[2]-2.0*f[3]+    f[4]),2)
                +coeff_weno5_gamma1*std::pow((3.0*f[2]-4.0*f[3]+    f[4]),2);
                
                
                q30= coeff_weno5_alpha[0][0]*f[0]
                +coeff_weno5_alpha[0][1]*f[1]
                +coeff_weno5_alpha[0][2]*f[2];
                
                q31= coeff_weno5_alpha[1][0]*f[1]
                +coeff_weno5_alpha[1][1]*f[2]
                +coeff_weno5_alpha[1][2]*f[3];
                
                q32= coeff_weno5_alpha[2][0]*f[2]
                +coeff_weno5_alpha[2][1]*f[3]
                +coeff_weno5_alpha[2][2]*f[4];
                
                alpha0=coeff_weno5_c[0]/std::pow((ss+is0),p);
                alpha1=coeff_weno5_c[1]/std::pow((ss+is1),p);
                alpha2=coeff_weno5_c[2]/std::pow((ss+is2),p);
                
                
                alphaSum=alpha0+alpha1+alpha2;
                
                omega0=alpha0/alphaSum;
                omega1=alpha1/alphaSum;
                omega2=alpha2/alphaSum;
                
                pOut[idx_rf]=omega0*q30+omega1*q31+omega2*q32;
            }
        }
    }
}

void nuc3d::weno5js::interpolationBoundaryR(const Field & fieldIN,
                                            const Field & boundaryR,
                                            const int uw,//uw=(1,-1)
                                            Field & fieldOUT,
                                            const int tilesize)
{
    switch(uw)
    {
        case 1:
            weno5jspBR(fieldIN,boundaryR, fieldOUT, tilesize);
            break;
        case -1:
            weno5jsnBR(fieldIN,boundaryR, fieldOUT, tilesize);
            break;
        default:
            std::cout<<"weno error: no such direction"<<std::endl;
            exit(-1);
    }
    
}


void nuc3d::weno5js::weno5jspBR(const Field & fieldIN,
                                const Field & boundaryR,
                                Field & fieldOUT,
                                const int tilesize)
{
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    double *pBND=boundaryR.getDataPtr();
    
    double omega0,omega1,omega2;
    double alpha0,alpha1,alpha2,alphaSum;
    double is0,is1,is2;
    double q30,q31,q32;
    
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
    
    double f[5];
    
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i+1;
                for(int z=-2;z<=2;z++)
                {
                    if((i+z)>=nx)
                    {
                        int idx_BND=nxBND*nyBND*k+nxBND*j+i+z-nx;
                        
                        f[2+z]=pBND[idx_BND];
                    }
                    else
                    {
                        int idx_f=nx*ny*k+nx*j+i+z;
                        f[2+z]=pIn[idx_f];
                    }
                }
                
                
                is0= coeff_weno5_gamma0*std::pow((    f[0]-2.0*f[1]+    f[2]),2)
                +coeff_weno5_gamma1*std::pow((    f[0]-4.0*f[1]+3.0*f[2]),2);
                
                is1= coeff_weno5_gamma0*std::pow((    f[1]-2.0*f[2]+    f[3]),2)
                +coeff_weno5_gamma1*std::pow((    f[1]-             f[3]),2);
                
                is2= coeff_weno5_gamma0*std::pow((    f[2]-2.0*f[3]+    f[4]),2)
                +coeff_weno5_gamma1*std::pow((3.0*f[2]-4.0*f[3]+    f[4]),2);
                
                
                q30= coeff_weno5_alpha[0][0]*f[0]
                +coeff_weno5_alpha[0][1]*f[1]
                +coeff_weno5_alpha[0][2]*f[2];
                
                q31= coeff_weno5_alpha[1][0]*f[1]
                +coeff_weno5_alpha[1][1]*f[2]
                +coeff_weno5_alpha[1][2]*f[3];
                
                q32= coeff_weno5_alpha[2][0]*f[2]
                +coeff_weno5_alpha[2][1]*f[3]
                +coeff_weno5_alpha[2][2]*f[4];
                
                alpha0=coeff_weno5_c[0]/std::pow((ss+is0),p);
                alpha1=coeff_weno5_c[1]/std::pow((ss+is1),p);
                alpha2=coeff_weno5_c[2]/std::pow((ss+is2),p);
                
                
                alphaSum=alpha0+alpha1+alpha2;
                
                omega0=alpha0/alphaSum;
                omega1=alpha1/alphaSum;
                omega2=alpha2/alphaSum;
                
                pOut[idx_rf]=omega0*q30+omega1*q31+omega2*q32;
            }
        }
    }
}

void nuc3d::weno5js::weno5jsnBR(const Field & fieldIN,
                                const Field & boundaryR,
                                Field & fieldOUT,
                                const int tilesize)
{
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    double *pBND=boundaryR.getDataPtr();
    
    double omega0,omega1,omega2;
    double alpha0,alpha1,alpha2,alphaSum;
    double is0,is1,is2;
    double q30,q31,q32;
    
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
    
    double f[5];
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i+1;
                
                for(int z=-2;z<=2;z++)
                {
                    if((i-z+1)>=nx)
                    {
                        int idx_BND=nxBND*nyBND*k+nxBND*j+i-z+1-nx;
                        
                        f[2+z]=pBND[idx_BND];
                    }
                    else
                    {
                        int idx_f=nx*ny*k+nx*j+i-z+1;
                        f[2+z]=pIn[idx_f];
                    }
                }
                
                
                is0= coeff_weno5_gamma0*std::pow((    f[0]-2.0*f[1]+    f[2]),2)
                +coeff_weno5_gamma1*std::pow((    f[0]-4.0*f[1]+3.0*f[2]),2);
                
                is1= coeff_weno5_gamma0*std::pow((    f[1]-2.0*f[2]+    f[3]),2)
                +coeff_weno5_gamma1*std::pow((    f[1]-             f[3]),2);
                
                is2= coeff_weno5_gamma0*std::pow((    f[2]-2.0*f[3]+    f[4]),2)
                +coeff_weno5_gamma1*std::pow((3.0*f[2]-4.0*f[3]+    f[4]),2);
                
                
                q30= coeff_weno5_alpha[0][0]*f[0]
                +coeff_weno5_alpha[0][1]*f[1]
                +coeff_weno5_alpha[0][2]*f[2];
                
                q31= coeff_weno5_alpha[1][0]*f[1]
                +coeff_weno5_alpha[1][1]*f[2]
                +coeff_weno5_alpha[1][2]*f[3];
                
                q32= coeff_weno5_alpha[2][0]*f[2]
                +coeff_weno5_alpha[2][1]*f[3]
                +coeff_weno5_alpha[2][2]*f[4];
                
                alpha0=coeff_weno5_c[0]/std::pow((ss+is0),p);
                alpha1=coeff_weno5_c[1]/std::pow((ss+is1),p);
                alpha2=coeff_weno5_c[2]/std::pow((ss+is2),p);
                
                
                alphaSum=alpha0+alpha1+alpha2;
                
                omega0=alpha0/alphaSum;
                omega1=alpha1/alphaSum;
                omega2=alpha2/alphaSum;
                
                pOut[idx_rf]=omega0*q30+omega1*q31+omega2*q32;
            }
        }
    }
}

#endif
