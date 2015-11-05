//
//  boundaryConditions.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/11/5.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#include "boundaryConditions.hpp"
#include "euler3d.h"


nuc3d::BufferMetric::BufferMetric(int nx,int ny,int nz):
jacobian(nx,ny,nz),
xi_xyz(3,Field(nx,ny,nz)),
eta_xyz(3,Field(nx,ny,nz)),
zeta_xyz(3,Field(nx,ny,nz))
{
    
}


nuc3d::BufferMetric::~BufferMetric()
{
    
}

nuc3d::boundaryCondition::boundaryCondition()
{
    
}

nuc3d::boundaryCondition::~boundaryCondition()
{
    
}

void nuc3d::boundaryCondition::initialBC(VectorBuffer &)
{
    
}


void nuc3d::boundaryCondition::setBC(VectorBuffer &,EulerData3D &)
{
    
}
void nuc3d::boundaryCondition::applyBC(VectorBuffer &,EulerData3D &)
{
    
}