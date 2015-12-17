#ifndef centraldifference6th_cpp
#define centraldifference6th_cpp
#include "centraldifference2nd.hpp"
#include "schemes.hpp"

nuc3d::centraldifference2nd::centraldifference2nd()
{
    
}

nuc3d::centraldifference2nd::~centraldifference2nd()
{
    
}

void nuc3d::centraldifference2nd::differentialInner(const Field & fieldIN,
                                                    const int dim0,
                                                    const int dim1,
                                                    const int dim2,
                                                    Field & fieldOUT,
                                                    const int tilesize)
{
    double flux_l;
    double flux_r;
    double h;
    
    int nx;
    int ny;
    int nz;
    
    nx=fieldIN.getSizeX();
    ny=fieldIN.getSizeY();
    nz=fieldIN.getSizeZ();
    
    int ibeg=tilesize*dim0;
    int iend=nx-tilesize*dim0;
    int jbeg=tilesize*dim1;
    int jend=ny-tilesize*dim1;
    int kbeg=tilesize*dim2;
    int kend=nz-tilesize*dim2;
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                
                int itempl=i-dim0;
                int jtempl=j-dim1;
                int ktempl=k-dim2;
                
                flux_l=fieldIN.getValue(itempl,jtempl,ktempl);
                
                int itempr=i+dim0;
                int jtempr=j+dim1;
                int ktempr=k+dim2;
                
                flux_r=fieldIN.getValue(itempr,jtempr,ktempr);
                
                h= nuc3d::cd2differential(flux_l,flux_r);
                
                fieldOUT.setValue(i,j,k,h);
            }
        }
    }
}


void nuc3d::centraldifference2nd::differentialBoundaryL(const Field & fieldIN,
                                                        const Field & boundaryL,
                                                        const int dim0,
                                                        const int dim1,
                                                        const int dim2,
                                                        Field & fieldOUT,
                                                        const int tilesize)
{
    double flux_l;
    double flux_r;
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
    
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                
                int itemp=i-dim0;
                int jtemp=j-dim1;
                int ktemp=k-dim2;
                
                if(itemp<0||jtemp<0||ktemp<0)
                {
                    int itemp1=itemp+dim0*boundaryL.getSizeX();
                    int jtemp1=jtemp+dim1*boundaryL.getSizeY();
                    int ktemp1=ktemp+dim2*boundaryL.getSizeZ();
                    
                    flux_l=boundaryL.getValue(itemp1,jtemp1,ktemp1);
                }
                else
                {
                    flux_l=fieldIN.getValue(itemp,jtemp,ktemp);
                }
                
                
                itemp=i+dim0;
                jtemp=j+dim1;
                ktemp=k+dim2;
                
                flux_r=fieldIN.getValue(itemp,jtemp,ktemp);
                
                h= nuc3d::cd2differential(flux_l,flux_r);
                fieldOUT.setValue(i,j,k,h);
            }
        }
    }
}


void nuc3d::centraldifference2nd::differentialBoundaryR(const Field & fieldIN,
                                                        const Field & boundaryR,
                                                        const int dim0,
                                                        const int dim1,
                                                        const int dim2,
                                                        Field & fieldOUT,
                                                        const int tilesize)
{
    double flux_l;
    double flux_r;
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
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int itemp=i-dim0;
                int jtemp=j-dim1;
                int ktemp=k-dim2;
                
                flux_l=fieldIN.getValue(itemp,jtemp,ktemp);
                
                
                itemp=i+dim0;
                jtemp=j+dim1;
                ktemp=k+dim2;
                
                if(itemp>=nx||jtemp>=ny||ktemp>=nz)
                {
                    int itemp1=itemp-dim0*nx;
                    int jtemp1=jtemp-dim1*ny;
                    int ktemp1=ktemp-dim2*nz;
                    
                    flux_r=boundaryR.getValue(itemp1,jtemp1,ktemp1);
                }
                else
                {
                    flux_r=fieldIN.getValue(itemp,jtemp,ktemp);
                }
                
                h= nuc3d::cd2differential(flux_l,flux_r);
                
                fieldOUT.setValue(i,j,k,h);
            }
        }
    }
}




#endif
