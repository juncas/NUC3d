#ifndef centraldifference2nd_cpp
#define centraldifference2nd_cpp
#include "centraldifference2nd.hpp"
#include "schemes.hpp"

nuc3d::centraldifference2nd::centraldifference2nd()
{
    
}

nuc3d::centraldifference2nd::~centraldifference2nd()
{
    
}

void nuc3d::centraldifference2nd::differentialInner(const Field & fieldIN,
                                                    Field & fieldOUT,
                                                    const int tilesize)
{
    
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    
    int nx=fieldIN.getSizeX();
    int ny=fieldIN.getSizeY();
    int nz=fieldIN.getSizeZ();
    
    int ibeg=tilesize;
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
                
                pOut[idx_f]=0.5*(pIn[idx_f+1]-pIn[idx_f-1]);
            }
        }
    }
}


void nuc3d::centraldifference2nd::differentialBoundaryL(const Field & fieldIN,
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
    
    int nxBND=boundaryL.getSizeX();
    int nyBND=boundaryL.getSizeY();
    int nzBND=boundaryL.getSizeZ();
    
    int ibeg=0;
    int iend=tilesize;
    int jbeg=0;
    int jend=ny;
    int kbeg=0;
    int kend=nz;
    
    double f[2];
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_f_out=nx*ny*k+nx*j+i;
                int idx_f=nx*ny*k+nx*j+i;
                int idx_BND=nxBND*nyBND*k+nxBND*j+(nxBND-1);
                
                f[0]=((i-1)<0)?pBND[idx_BND]:pIn[idx_f-1];
                f[1]=pIn[idx_f+1];
                
                pOut[idx_f_out]=0.5*(f[1]-f[0]);
                
            }
        }
    }
}


void nuc3d::centraldifference2nd::differentialBoundaryR(const Field & fieldIN,
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
    
    int nxBND=boundaryR.getSizeX();
    int nyBND=boundaryR.getSizeY();
    int nzBND=boundaryR.getSizeZ();
    
    int ibeg=nx-tilesize;
    int iend=nx;
    int jbeg=0;
    int jend=ny;
    int kbeg=0;
    int kend=nz;
    
    double f[2];
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_f_out=nx*ny*k+nx*j+i;
                int idx_f=nx*ny*k+nx*j+i;
                int idx_BND=nxBND*nyBND*k+nxBND*j;
                
                f[0]=pIn[idx_f-1];
                f[1]=((i+1)>=nx)?pBND[idx_BND]:pIn[idx_f+1];
                
                pOut[idx_f_out]=0.5*(f[1]-f[0]);
            }
        }
    }
}


#endif
