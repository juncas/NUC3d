//
//  singleBlock.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/10/20.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#include "singleBlock.h"

nuc3d::singleBlock::singleBlock(int& argc, char **&argv):
block(),
myComm(argc, argv),
myCtrler(),
myPhys(),
myOperator()
{
    initialBlock();
}

nuc3d::singleBlock::~singleBlock()
{
    
}

void nuc3d::singleBlock::initialBlock()
{
    std::stringstream ss;
    int ProcId = myComm.getMyId();
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
    
    myFile.open(filename_mesh);
    
    if (myFile)
    {
        myFile >> nx0 >> ny0 >> nz0;
        
        block::initial(nx0, ny0, nz0, myPhys);
        myOperator.initial(myFluxes->getPrimatives());
        
        double x,y,z;
        
        for(int k=0;k<nz+1;k++)
        {
            for(int j=0;j<ny+1;j++)
            {
                for(int i=0;i<nx+1;i++)
                {
                    myFile>>x>>y>>z;
                    xyz[0].setValue(i, j, k, x);
                    xyz[1].setValue(i, j, k, y);
                    xyz[2].setValue(i, j, k, z);
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
    
    myFile.open(filename_flow);
    
    if (myFile)
    {
        myFile >> nx0 >> ny0 >> nz0;
        if (nx0==nx&&ny0==ny&&nz0==nz)
        {
            double value;
            for (auto iter=(myPDE.getQ()).begin();
                 iter!=(myPDE.getQ()).end(); iter++)
            {
                for(int k=0;k<nz+1;k++)
                {
                    for(int j=0;j<ny+1;j++)
                    {
                        for(int i=0;i<nx+1;i++)
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

void nuc3d::singleBlock::initial()
{
    myPhys.initial(myPDE,myFluxes);
}

void nuc3d::singleBlock::solve()
{
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
            solveIntegral();
        }
        
        printRES();
        
        if (myCtrler.ifsave())
        {
            output();
        }

        myCtrler.increaseStep();
    }
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
