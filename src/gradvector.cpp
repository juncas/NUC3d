//
//  gradvector.cpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/12.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//

#include "gradvector.hpp"

nuc3d::gradvector::gradvector(int nx0,int ny0,int nz0):
f_xi(nx0,ny0,nz0),
f_eta(ny0,nz0,nx0),
f_zeta(nz0,nx0,ny0),
dxi(nx0,ny0,nz0),
deta(ny0,nz0,nx0),
dzeta(nz0,nx0,ny0)
{
    
}

void nuc3d::gradvector::setGrad(Field &myField)
{
    double *pField=myField.getDataPtr();
    double *pf_xi=f_xi.getDataPtr();
    double *pf_eta=f_eta.getDataPtr();
    double *pf_zeta=f_zeta.getDataPtr();
    
    int nx0=myField.getSizeX();
    int ny0=myField.getSizeY();
    int nz0=myField.getSizeZ();
    
    for (int k=0; k<nz0; k++)
    {
        for (int j=0; j<ny0; j++)
        {
            for (int i=0; i<nx0; i++)
            {
                int idx=nx0*ny0*k+nx0*j+i;
                pf_xi[idx]=pField[idx];
            }
        }
    }
    
    for (int k=0; k<nz0; k++)
    {
        for (int j=0; j<ny0; j++)
        {
            for (int i=0; i<nx0; i++)
            {
                int idx=nx0*ny0*k+nx0*j+i;
                int idx_eta=ny0*nz0*i+ny0*k+j;
                
                pf_eta[idx_eta]=pField[idx];
            }
        }
    }
    
    for (int k=0; k<nz0; k++)
    {
        for (int j=0; j<ny0; j++)
        {
            for (int i=0; i<nx0; i++)
            {
                int idx=nx0*ny0*k+nx0*j+i;
                int idx_zeta=nz0*nx0*j+nz0*i+k;
                
                pf_zeta[idx_zeta]=pField[idx];
            }
        }
    }
    
}
nuc3d::gradvector::~gradvector()
{}
