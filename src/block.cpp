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

void nuc3d::block::initial(fieldOperator3d &myOP,
		physicsModel &myPhyMod,
		MPIComunicator3d_nonblocking &myMPI,
		boundaryCondition &myBC,
		IOController &myIO)

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

	myFile.open(filename_mesh);
	if (myFile)
	{
		myFile >> nx0 >> ny0 >> nz0>> bfsize;
		initialData(nx0, ny0, nz0, myPhyMod);
		myBC.initialBC(mybuffer,myMPI);
	}
	myFile.close();


}

void nuc3d::block::initialData(int nx0,int ny0,int nz0,physicsModel &myPhy)
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
		//        myFluxes=std::make_shared<EulerReactiveData3D>(nx0,ny0,nz0,myPhy.getEqNum());
	}
	else if ("NaiverStokes3d"==myPhy.getMyModelName())
	{
		//        myFluxes=std::make_shared<NaiverStokesData3d>(nx0,ny0,nz0,myPhy.getEqNum());
	}
	else if ("NaiverStokesReactive3d"==myPhy.getMyModelName())
	{
		//        myFluxes=std::make_shared<NaiverStokesReactiveData3d>(nx0,ny0,nz0,myPhy.getEqNum());
	}
	else
	{
		std::cout <<"Model Name:"<<myPhy.getMyModelName()
			<<"does not exist!"
			<<std::endl;
		exit(-1);
	}
	
	for(int n=0;n<myPhy.getEqNum();n++)
		mybuffer.push_back(bufferData(nx,ny,nz,bfsize));

	std::cout<<"Field data has been allocated!"<<std::endl;
}

void nuc3d::block::solve(fieldOperator3d &myOP,
		physicsModel &myPhyMod,
		MPIComunicator3d_nonblocking &myMPI,
		boundaryCondition &myBC,
		IOController &myIO)
{
	int step=0;

	while (myOP.getSteps()!=step)
	{
		myFluxes->solve(myPDE, myOP, mybuffer, myPhyMod, myMPI, myBC);
		myPDE.solve(myOP, step++);
	}

	istep++;
	myIO.myIOController["currentStep"]=istep;

	dt=myPDE.getDt();
	time+=dt;

	myIO.myTimeController["dt"]=dt;
	myIO.myTimeController["currentTime"]=time;


	RES=myPDE.getRes();
}

void nuc3d::block::printStatus()
{
	std::cout<<"step "<<istep<< ", time = "<<time<<", residual = "<<RES<<std::endl;
}

void nuc3d::block::ReadXYZ()
{

	myFluxes->initialxyz(xyz);
}

void nuc3d::block::ReadPDE()
{

}

