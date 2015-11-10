//
//  PDEData3d.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/11/4.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#include "PDEData3d.hpp"
#include "mpi.h"
/**************************************************************************************
 Member functions of class: PDEData3d
 **************************************************************************************/
nuc3d::PDEData3d::PDEData3d()
{
    
}

nuc3d::PDEData3d::PDEData3d(int nx0,int ny0,int nz0,int neqs):
nEquations(neqs),
Q(neqs,Field(nx0,ny0,nz0)),
Q_Euler(neqs,Field(nx0,ny0,nz0)),
Q_work(neqs,Field(nx0,ny0,nz0)),
dt_local(0.0),
dt_global(0.0),
res_global(0.0),
res_local(0.0)
{
    
}

void nuc3d::PDEData3d::initPDEData3d(int nx0,int ny0,int nz0,int neqs)
{
    nEquations=neqs;
    for(int i=0;i<neqs;i++)
    {
        Q_Euler.push_back(Field(nx0,ny0,nz0));
        Q_work.push_back(Field(nx0,ny0,nz0));
    }
    
    dt_local=0.0;
    dt_global=0.0;
    res_global=0.0;
    res_local=0.0;

    
}


nuc3d::VectorField& nuc3d::PDEData3d::getRHS()
{
    return RHS;
}

nuc3d::VectorField& nuc3d::PDEData3d::getQ()
{
    return Q_work;
}

void nuc3d::PDEData3d::solve(fieldOperator3d &myOP,
                             int step)
{
    if(step==1)
    {
        Q_Euler=Q_work;
        
        MPI_Allreduce(&dt_local, &dt_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    }
    myOP.timeIntegral(RHS, Q_Euler,Q_work, dt_global, step);
    
    if(step==3)
    {
        setRES();
    }
}

void  nuc3d::PDEData3d::setRES()
{
    auto beg=Q_Euler.begin();
    auto end=Q_Euler.end();
    
    int nx=beg->getSizeX();
    int ny=beg->getSizeY();
    int nz=beg->getSizeZ();
    
    res_local=0.0;
    for(auto iter=beg;iter!=end;iter++)
    {
        for(int k=0;k<nz;++k)
        {
            for(int j=0;j<ny;++j)
            {
                for(int i=0;i<nx;++i)
                {
                    double u=Q_work[iter-beg].getValue(i, j, k);
                    double u0=iter->getValue(i, j, k);
                    res_local+=std::pow(u-u0,2);
                }
            }
        }
    }
    
    res_local=std::sqrt(res_local/(nx*ny*nz*Q_work.size()));
    
    MPI_Allreduce(&res_local, &res_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

nuc3d::PDEData3d::~PDEData3d()
{}