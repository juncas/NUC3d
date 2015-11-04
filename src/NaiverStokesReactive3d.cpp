//
//  NaiverStokesReactive3d.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/11/3.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#include "NaiverStokesReactive3d.h"
/**************************************************************************************
 Member functions of class: NaiverStokesReactiveData3d
 **************************************************************************************/
nuc3d::NaiverStokesReactiveData3d::NaiverStokesReactiveData3d(int nx0, int ny0, int nz0, int neqs):
EulerData3D(nx0,ny0,nz0,neqs),
EulerReactiveData3D(nx0,ny0,nz0,neqs),
NaiverStokesData3d(nx0,ny0,nz0,neqs)
{
    
}

nuc3d::VectorField& nuc3d::NaiverStokesReactiveData3d::getDrivativeXi()
{
    for(auto iter=dfdxi.begin();iter!=dfdxi.end();iter++)
    {
        *iter-=dfvdxi[iter-dfdxi.begin()];
        *iter-source_xi[iter-dfdxi.begin()];
    }
    
    return this->dfdxi;
}

nuc3d::VectorField& nuc3d::NaiverStokesReactiveData3d::getDrivativeEta()
{
    for(auto iter=dgdeta.begin();iter!=dgdeta.end();iter++)
    {
        *iter-=dgdeta[iter-dfdxi.begin()];
        *iter-=source_eta[iter-dfdxi.begin()];
    }
    
    return this->dgdeta;
}

nuc3d::VectorField& nuc3d::NaiverStokesReactiveData3d::getDrivativeZeta()
{
    for(auto iter=dhdzeta.begin();iter!=dhdzeta.end();iter++)
    {
        *iter-=dhvdzeta[iter-dfdxi.begin()];
        *iter-=source_zeta[iter-dfdxi.begin()];
    }
    
    return this->dhdzeta;
}

void nuc3d::NaiverStokesReactiveData3d::solveLocal()
{
    EulerReactiveData3D::solveLocal();
    NaiverStokesReactiveData3d::solveLocal();
}

nuc3d::NaiverStokesReactiveData3d::~NaiverStokesReactiveData3d()
{}

