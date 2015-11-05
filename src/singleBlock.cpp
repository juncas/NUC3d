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
myOperator()
{
    
}

nuc3d::singleBlock::~singleBlock()
{
    
}

/*==========================================================================

 Initial Block functions:
    -initialBlock()
 =========================================================================*/

void nuc3d::singleBlock::initialBlock()
{
    std::stringstream ss;
    int ProcId = myMPI.getMyId();
    int nx0, ny0, nz0;
    std::string forename_mesh = ("mesh_");
    std::string forename_flow = ("flow_");
    std::string midname;
    std::string tailname = (".dat");
    
    ss << ProcId;
    ss >> midname;
    
    std::string filename_mesh = forename_mesh + midname + tailname;
    std::string filename_flow = forename_flow + midname + tailname;
    
    std::ifstream myFile;
    
    
//read grid data
    myFile.open(filename_mesh);
    
    if (myFile)
    {
        myFile >> nx0 >> ny0 >> nz0;
        
    //initial block
        myBlock.initial(nx0, ny0, nz0, myPhys);
        
    //initial field operator
        myOperator.initial(myBlock.myFluxes->getPrimatives());
        
    //initial buffer
        for (int i=0; i<myPhys.getEqNum(); i++)
        {
            myBlock.mybuffer.push_back(bufferData(nx0, ny0, nz0, myOperator.getBufferSize()));
        }
        
        
    //read location data x,y,z at cell corner
        double x,y,z;
        
        for(int k=0;k<myBlock.nz+1;k++)
        {
            for(int j=0;j<myBlock.ny+1;j++)
            {
                for(int i=0;i<myBlock.nx+1;i++)
                {
                    myFile>>x>>y>>z;
                    myBlock.xyz[0].setValue(i, j, k, x);
                    myBlock.xyz[1].setValue(i, j, k, y);
                    myBlock.xyz[2].setValue(i, j, k, z);
                }
            }
        }
        
        
    }
    else
    {
        std::cout << "file " << filename_mesh << " does not exist" << std::endl;
        myComm.AbortMPI();
    }
    myFile.close();
    
    
//read flow data
    myFile.open(filename_flow);
    
    if (myFile)
    {
        myFile >> nx0 >> ny0 >> nz0;
        if (nx0==myBlock.nx&&ny0==myBlock.ny&&nz0==myBlock.nz)
        {
            double value;
            for (auto iter=(myBlock.myPDE.getQ()).begin();
                 iter!=(myBlock.myPDE.getQ()).end(); iter++)
            {
                for(int k=0;k<myBlock.nz+1;k++)
                {
                    for(int j=0;j<myBlock.ny+1;j++)
                    {
                        for(int i=0;i<myBlock.nx+1;i++)
                        {
                            myFile>>value;
                            iter->setValue(i, j, k, value);
                        }
                    }
                }
            }
            
        }
        else
        {
            std::cout<<"grid number does not match data number"
            <<std::endl;
        }
    }
    else
    {
        std::cout << "file " << filename_mesh << " does not exist" << std::endl;
        myComm.AbortMPI();
    }
    
    myFile.close();
}

/*==========================================================================
 
 Block solve functions:
    -solvePDE()
        --solveRiemann()
        --solve BoundaryConditions()
        --solveInvicidFlux()
        --solveViscousFlux()
        --solveGetRHS()
        --solveIntegral(nstep)
 =========================================================================*/

void nuc3d::singleBlock::solvePDE()
{
    initialPDE();
    while (myCtrler.ifsolve())
    {
        int nstep=0;
        while (nstep<myOperator.getSteps())
        {
            solveRiemann();
            solveBoundaryConditions();
            solveInvicidFlux();
            solveViscousFLux();
            solveGetRHS();
            solveIntegral(nstep);
            nstep++;
        }
        
        printRES();
        
        if (myCtrler.ifsave())
        {
            output();
        }

        myCtrler.renew();
    }
}

void nuc3d::singleBlock::initialPDE()
{
    myPhys.initial(myBlock.myPDE,myBlock.myFluxes);
}


void nuc3d::singleBlock::solveRiemann()
{
    myPhys.solve(myBlock.myPDE,myBlock.myFluxes);
}

void nuc3d::singleBlock::solveBoundaryConditions()
{
    
}

void nuc3d::singleBlock::output()
{
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
