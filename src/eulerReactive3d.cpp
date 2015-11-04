//
//  eulerReactive3d.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/11/3.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#include "eulerReactive3d.h"

/**************************************************************************************
 Member functions of class: EulerReactiveData3D
 **************************************************************************************/
nuc3d::EulerReactiveData3D::EulerReactiveData3D( int nx0, int ny0, int nz0, int neqs):
EulerData3D(nx0,ny0,nz0,neqs),
source_xi(neqs,Field(nx0,ny0,nz0,0.0)),
source_eta(neqs,Field(nx0,ny0,nz0,0.0)),
source_zeta(neqs,Field(nx0,ny0,nz0,0.0))
{
    
}


nuc3d::VectorField& nuc3d::EulerReactiveData3D::getDrivativeXi()
{
    for(auto iter=dfdxi.begin();iter!=dfdxi.end();iter++)
        *iter-=source_xi[iter-dfdxi.begin()];
    
    return this->dfdxi;
}

nuc3d::VectorField& nuc3d::EulerReactiveData3D::getDrivativeEta()
{
    for(auto iter=dgdeta.begin();iter!=dgdeta.end();iter++)
        *iter-=source_eta[iter-dfdxi.begin()];
    
    return this->dgdeta;
}

nuc3d::VectorField& nuc3d::EulerReactiveData3D::getDrivativeZeta()
{
    for(auto iter=dhdzeta.begin();iter!=dhdzeta.end();iter++)
        *iter-=source_zeta[iter-dfdxi.begin()];
    
    return this->dhdzeta;
}

void nuc3d::EulerReactiveData3D::solveLocal()
{
    setSource();
}

void nuc3d::EulerReactiveData3D::setSource()
{
    
}

nuc3d::EulerReactiveData3D::~EulerReactiveData3D()
{
    
}

