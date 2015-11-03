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
        for (int i=0; i<myPhys.getEqNum(); i++)
        {
            mybuffer.push_back(bufferData(nx0, ny0, nz0, myOperator.getBufferSize()));
        }
        
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

void nuc3d::singleBlock::initialPDE()
{
    myPhys.initial(myPDE,myFluxes);
}

void nuc3d::singleBlock::solvePDE()
{
    while (myCtrler.ifsolve())
    {
        int nstep=0;
        while (nstep<myOperator.getSteps())
        {
            solveRiemann();
            solveBoundaryConditions();
            solveInvicidFlux();
            solveViscousFLux(myFluxesMap[myFluxName]);
            solveGetRHS(myPDE);
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

void nuc3d::singleBlock::output()
{
}


void nuc3d::singleBlock::solveRiemann()
{
    myPhys.solve(block::myPDE, myFluxes);
}

void nuc3d::singleBlock::solveBoundaryConditions()
{
    
}

void nuc3d::singleBlock::solveInvicidFlux()
{
	solveInvicidFluxL(myFluxes->getFluxXi(), mybuffer, 0);
	solveInvicidFluxL(myFluxes->getFluxEta(), mybuffer, 1);
	solveInvicidFluxL(myFluxes->getFluxZeta(), mybuffer, 2);

	solveInvicidFluxR(myFluxes->getFluxXi(), mybuffer, 0);
	solveInvicidFluxR(myFluxes->getFluxEta(), mybuffer, 1);
	solveInvicidFluxR(myFluxes->getFluxZeta(), mybuffer, 2);
    
    myFluxes->setDerivativesInv();
}

void nuc3d::singleBlock::solveInvicidFluxL(EulerFlux &myFlux,std::vector<bufferData> &myBuff, int dir)
{
	VectorField &pFlux = myFlux.FluxL;
	VectorField &pReconFlux = myFlux.reconstFluxL;

	for (auto iter = pFlux.begin(); iter != pFlux.end(); iter++)
	{
		Field &rf = pReconFlux[iter - pFlux.begin()];
		bufferData &bf = myBuff[iter - pFlux.begin()];

		bf.setBufferSend(*iter);

		myComm.bufferSendRecv(bf, 0);

		myOperator.reconstructionInner(*iter, 0, 1, rf);

		myComm.waitAllSendRecv(bf);

		myOperator.reconstructionBoundary(*iter, bf.BufferRecv[dir], bf.BufferRecv[dir+1], dir, 1, rf);
	
		MPI_Barrier();
	}
    
}

void nuc3d::singleBlock::solveInvicidFluxR(EulerFlux &myFlux, std::vector<bufferData> &myBuff, int dir)
{
	VectorField &pFlux = myFlux.FluxR;
	VectorField &pReconFlux = myFlux.reconstFluxR;

	for (auto iter = pFlux.begin(); iter != pFlux.end(); iter++)
	{
		Field &rf = pReconFlux[iter - pFlux.begin()];
		bufferData &bf = myBuff[iter - pFlux.begin()];

		bf.setBufferSend(*iter);

		myComm.bufferSendRecv(bf, 0);

		myOperator.reconstructionInner(*iter, 0, 1, rf);

		myComm.waitAllSendRecv(bf);

		myOperator.reconstructionBoundary(*iter, bf.BufferRecv[dir], bf.BufferRecv[dir + 1], dir, -1, rf);

		MPI_Barrier();
	}

}

void nuc3d::singleBlock::solveViscousFLux(EulerData3D &myEuler  )
{
    
}

void nuc3d::singleBlock::solveViscousFLux(EulerReactiveData3D &myEuler  )
{
    
}


void nuc3d::singleBlock::solveViscousFLux(NaiverStokesData3d &myEuler  )
{
    myEuler.du
}


void nuc3d::singleBlock::solveViscousFLux(NaiverStokesReactiveData3d &myEuler  )
{
}

void nuc3d::singleBlock::solveGetRHS()
{
    myPDE.getRHS(*myFluxes);
}

void nuc3d::singleBlock::solveIntegral(int step)
{
	myOperator.timeIntegral(myPDE.getRHS(*myFluxes),
                            myPDE.getQ(),
                            myCtrler.getValue("dt"),
                            step);
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
