//
//  block.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/10/20.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#include "block.hpp"

nuc3d::block::block(int nx0,int ny0,int nz0,physicsModel &myPhy):
nx(nx0),
ny(ny0),
nz(nz0),
xyz(3,Field(nx0+1,ny0+1,nz0+1)),
xyz_center(3,Field(nx0,ny0,nz0)),
myPDE(nx0,ny0,nz0,myPhy.getEqNum())
{
    if("Euler3d"==myPhy.getMyModelName())
    {
        myFluxes=std::make_shared<nuc3d::EulerData3D>(nx0,ny0,nz0,myPhy.getEqNum());
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

nuc3d::block::~block()
{
    
}