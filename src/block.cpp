//
//  block.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/10/20.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#include "block.h"
#include "euler3d.h"
#include "eulerReactive3d.h"
#include "NaiverStokes3d.h"
#include "NaiverStokesReactive3d.h"
#include "physicsModel.h"
#include "MPICommunicator.h"
#include "IOcontroller.h"
#include "boundaryConditions.hpp"

 nuc3d::block::block()
{}

nuc3d::block::~block()
{}

void nuc3d::block::initial(int nx0,int ny0,int nz0,physicsModel &myPhy)
{
    nx=nx0;
    ny=ny0;
    nz=nz0;
    
    for(int i=0;i<3;i++)
    {
        xyz.push_back(Field(nx0+1,ny0+1,nz0+1));
        xyz_center.push_back(Field(nx0,ny0,nz0));
    }
    
    myPDE.initPDEData3d(nx, ny, nz, myPhy.getEqNum());
    
    if("Euler3d"==myPhy.getMyModelName())
    {
        myFluxes=std::make_shared<EulerData3D>(nx0,ny0,nz0,myPhy.getEqNum());
    }
    else if ("EulerReactive3d"==myPhy.getMyModelName())
    {
        myFluxes=std::make_shared<EulerReactiveData3D>(nx0,ny0,nz0,myPhy.getEqNum());
    }
    else if ("NaiverStokes3d"==myPhy.getMyModelName())
    {
        myFluxes=std::make_shared<NaiverStokesData3d>(nx0,ny0,nz0,myPhy.getEqNum());
    }
    else if ("NaiverStokesReactive3d"==myPhy.getMyModelName())
    {
        myFluxes=std::make_shared<NaiverStokesReactiveData3d>(nx0,ny0,nz0,myPhy.getEqNum());
    }
    else
    {
        std::cout <<"Model Name:"<<myPhy.getMyModelName()
        <<"does not exist!"
        <<std::endl;
        exit(-1);
    }
    
}

void nuc3d::block::solve(fieldOperator3d &myOP,
                         physicsModel &myPhyMod,
                         MPIComunicator3d_nonblocking &myMPI,
                         boundaryCondition &myBC,
                         int step)
{
    myFluxes->solve(myPDE, myOP, mybuffer, myPhyMod, myMPI,myBC);
    myPDE.solve(myOP, step);
}
