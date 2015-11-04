//
//  PDEData3d.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/11/4.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#include "PDEData3d.hpp"
/**************************************************************************************
 Member functions of class: PDEData3d
 **************************************************************************************/
nuc3d::PDEData3d::PDEData3d()
{
    
}

nuc3d::PDEData3d::PDEData3d(int nx0,int ny0,int nz0,int neqs):
nEquations(neqs),
Q_Euler(neqs,Field(nx0,ny0,nz0))
{
    
}

void nuc3d::PDEData3d::initPDEData3d(int nx0,int ny0,int nz0,int neqs)
{
    nEquations=neqs;
    for(int i=0;i<neqs;i++)
    {
        Q_Euler.push_back(Field(nx0,ny0,nz0));
    }
    
}


nuc3d::VectorField& nuc3d::PDEData3d::getRHS()
{
    return RHS;
}

nuc3d::VectorField& nuc3d::PDEData3d::getQ()
{
    return Q_Euler;
}

nuc3d::PDEData3d::~PDEData3d()
{}