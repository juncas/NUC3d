//
//  singleBlock.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/10/20.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//
#include "singleBlock.h"

nuc3d::singleBlock::singleBlock():
myBlock(),
myMPI(),
myIO(),
myPhys(),
myOperator(),
myBC()
{
    
}

nuc3d::singleBlock::~singleBlock()
{
    
}

/*==========================================================================
 
 Initial Block functions:
 -initialBlock()
 =========================================================================*/

void nuc3d::singleBlock::loop()
{
    initialBlock();
    MPI_Barrier(MPI_COMM_WORLD);
    solvePDE();
}

void nuc3d::singleBlock::initialBlock()
{
    if(myMPI.getMyId()==0) std::cout<<"Initializing......\n";
    myBlock.initial(myOperator,myPhys,myMPI,myBC,myIO);
}

/*==========================================================================
 
 Block solve functions:
 -solvePDE()
 =========================================================================*/

void nuc3d::singleBlock::solvePDE()
{
    if(myMPI.getMyId()==0) std::cout<<"Start solving......\n";
    while (myIO.ifsolve(myBlock.getTime(),myBlock.getStep()))
    {
        myBlock.solve(myOperator, myPhys, myMPI, myBC, myIO);
        
        if(0==myMPI.getMyId()) myBlock.printStatus();
        
        if(myIO.ifpost(myBlock.getStep()))
        {
            postprocess();
        }
        
        if(myIO.ifsave(myBlock.getStep()))
        {
            output();
        }
    }
    
    postprocess();
    output();
    
    if(0==myMPI.getMyId()) std::cout<<"Calculation finished!!"<<std::endl;
}

void nuc3d::singleBlock::postprocess()
{
    if(0==myMPI.getMyId()) std::cout<<"Post Processing..."<<std::endl;    myBlock.Post(myOperator, myPhys, myMPI, myBC, myIO);
    myIO.renewIOcontroller();
}

void nuc3d::singleBlock::output()
{
    if(0==myMPI.getMyId()) std::cout<<"Saving..."<<std::endl;
    myBlock.Output(myOperator, myPhys, myMPI, myBC, myIO);
    myIO.renewIOcontroller();
}


void nuc3d::singleBlock::readData(std::ifstream &myFile, VectorField &dataVec)
{
    double value;
    
    int nx=(dataVec.begin())->getSizeX();
    int ny=(dataVec.begin())->getSizeY();
    int nz=(dataVec.begin())->getSizeZ();
    
    for (auto iter = dataVec.begin(); iter != dataVec.end(); ++iter)
    {
        for(int k=0;k<nz;k++)
        {
            for(int j=0;j<ny;j++)
            {
                for(int i=0;i<nx;i++)
                {
                    myFile>>value;
                    iter->setValue(i, j, k, value);
                }
            }
        }
    }
}

void nuc3d::singleBlock::writeData(std::ofstream &myFile, VectorField &dataVec)
{
    double value;
    int nx=(dataVec.begin())->getSizeX();
    int ny=(dataVec.begin())->getSizeY();
    int nz=(dataVec.begin())->getSizeZ();
    
    for (auto iter = dataVec.begin(); iter != dataVec.end(); ++iter)
    {
        for(int k=0;k<nz;k++)
        {
            for(int j=0;j<ny;j++)
            {
                for(int i=0;i<nx;i++)
                {
                    value=iter->getValue(i,j,k);
                    myFile<<value<<std::endl;
                }
            }
        }
    }
    
}
