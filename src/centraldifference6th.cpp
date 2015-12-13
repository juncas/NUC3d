#ifndef centraldifference6th_cpp
#define centraldifference6th_cpp
#include "centraldifference6th.h"

nuc3d::centraldifference6th::centraldifference6th()
{
    
}

nuc3d::centraldifference6th::~centraldifference6th()
{
    
}

double nuc3d::centraldifference6th::cd6differential(const double *fl,
                                                    const double *fr)
{
    double df;
    
    df = coeff_cd6_alpha[2]*(fr[2]-fl[2])
    +coeff_cd6_alpha[1]*(fr[1]-fl[1])
    +coeff_cd6_alpha[0]*(fr[0]-fl[0]);
    
    return df;
}

void nuc3d::centraldifference6th::differentialInner(const Field & fieldIN,
                                                    const int dim0,
                                                    const int dim1,
                                                    const int dim2,
                                                    Field & fieldOUT)
{
    double flux_l[3];
    double flux_r[3];
    double h;
    
    int nx;
    int ny;
    int nz;
    
    nx=fieldIN.getSizeX();
    ny=fieldIN.getSizeY();
    nz=fieldIN.getSizeZ();
    
    int strideL=3;
    int strideR=3;
    
    int ibeg=strideL*dim0;
    int iend=nx-strideR*dim0;
    int jbeg=strideL*dim1;
    int jend=ny-strideR*dim1;
    int kbeg=strideL*dim2;
    int kend=nz-strideR*dim2;
    
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
                    for(int z=0;z<3;z++)
                    {
                        int stride=-(z+1);
                        int itemp=i+stride*dim0;
                        int jtemp=j+stride*dim1;
                        int ktemp=k+stride*dim2;
                        
                        flux_l[z]=fieldIN.getValue(itemp,jtemp,ktemp);
                    }
                    
                    for(int z=0;z<3;z++)
                    {
                        int stride=z+1;
                        int itemp=i+stride*dim0;
                        int jtemp=j+stride*dim1;
                        int ktemp=k+stride*dim2;
                        
                        flux_r[z]=fieldIN.getValue(itemp,jtemp,ktemp);
                    }
                    
                    h= cd6differential(flux_l,flux_r);
                    
                    fieldOUT.setValue(i,j,k,h);
                }
            }
        }
    }
}


void nuc3d::centraldifference6th::differentialBoundaryL(const Field & fieldIN,
                                                        const Field & boundaryL,
                                                        const int dim0,
                                                        const int dim1,
                                                        const int dim2,
                                                        Field & fieldOUT
                                                        )
{
    double flux_l[3];
    double flux_r[3];
    double h;
    
    int nx;
    int ny;
    int nz;
    
    nx=fieldIN.getSizeX();
    ny=fieldIN.getSizeY();
    nz=fieldIN.getSizeZ();
    
    int strideL=3;
    
    int ibeg=0;
    int iend=(strideL-nx)*dim0+nx;
    int jbeg=0;
    int jend=(strideL-ny)*dim1+ny;
    int kbeg=0;
    int kend=(strideL-nz)*dim2+nz;
    
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
                    for(int z=0;z<3;z++)
                    {
                        int stride=-(z+1);
                        int itemp=i+stride*dim0;
                        int jtemp=j+stride*dim1;
                        int ktemp=k+stride*dim2;
                        
                        if(itemp<0||jtemp<0||ktemp<0)
                        {
                            int itemp1=itemp+dim0*boundaryL.getSizeX();
                            int jtemp1=jtemp+dim1*boundaryL.getSizeY();
                            int ktemp1=ktemp+dim2*boundaryL.getSizeZ();
                            
                            flux_l[z]=boundaryL.getValue(itemp1,jtemp1,ktemp1);
                        }
                        else
                        {
                            flux_l[z]=fieldIN.getValue(itemp,jtemp,ktemp);
                        }
                    }
                    
                    for(int z=0;z<3;z++)
                    {
                        int stride=z+1;
                        int itemp=i+stride*dim0;
                        int jtemp=j+stride*dim1;
                        int ktemp=k+stride*dim2;
                        
                        flux_r[z]=fieldIN.getValue(itemp,jtemp,ktemp);
                    }
                    
                    h= cd6differential(flux_l,flux_r);
                    
                    fieldOUT.setValue(i,j,k,h);
                }
            }
        }
    }
}


void nuc3d::centraldifference6th::differentialBoundaryR(const Field & fieldIN,
                                                        const Field & boundaryR,
                                                        const int dim0,
                                                        const int dim1,
                                                        const int dim2,
                                                        Field & fieldOUT
                                                        )
{
    double flux_l[3];
    double flux_r[3];
    double h;
    
    int nx;
    int ny;
    int nz;
    
    nx=fieldIN.getSizeX();
    ny=fieldIN.getSizeY();
    nz=fieldIN.getSizeZ();
    
    int strideR=3;
    
    int ibeg=(nx-strideR)*dim0;
    int iend=nx;
    int jbeg=(ny-strideR)*dim1;
    int jend=ny;
    int kbeg=(nz-strideR)*dim2;
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
                    for(int z=0;z<3;z++)
                    {
                        int stride=-(z+1);
                        int itemp=i+stride*dim0;
                        int jtemp=j+stride*dim1;
                        int ktemp=k+stride*dim2;
                        
                        flux_l[z]=fieldIN.getValue(itemp,jtemp,ktemp);
                    }
                    
                    for(int z=0;z<3;z++)
                    {
                        int stride=z+1;
                        int itemp=i+stride*dim0;
                        int jtemp=j+stride*dim1;
                        int ktemp=k+stride*dim2;
                        
                        if(itemp>=nx||jtemp>=ny||ktemp>=nz)
                        {
                            int itemp1=itemp-dim0*nx;
                            int jtemp1=jtemp-dim1*ny;
                            int ktemp1=ktemp-dim2*nz;
                            
                            flux_r[z]=boundaryR.getValue(itemp1,jtemp1,ktemp1);
                        }
                        else
                        {
                            flux_r[z]=fieldIN.getValue(itemp,jtemp,ktemp);
                        }

                    }

                    
                    h= cd6differential(flux_l,flux_r);
                    
                    fieldOUT.setValue(i,j,k,h);
                }
            }
        }
    }
}


#endif
