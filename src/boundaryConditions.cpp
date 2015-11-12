//
//  boundaryConditions.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/11/5.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#include "boundaryConditions.hpp"
#include "MPICommunicator.h"
#include "euler3d.h"
#include "PDEData3d.hpp"
#include "physicsModel.h"

nuc3d::boundaryCondition::boundaryCondition()
{
}

nuc3d::boundaryCondition::~boundaryCondition()
{

}

const static std::string BCTypes[2]={"Exterior","Interior"};
const static std::string BCnames[5]={"Inlet condition",
	"Outlet condition",
	"Wall condition",
	"Symmetric condition",
	"Periodic condition"};
void nuc3d::boundaryCondition::initialBC(VectorBuffer &myBuffer,
		MPIComunicator3d_nonblocking &myMPI)
{
	std::stringstream ss;
	int ProcId = myMPI.getMyId();
	std::string forename_bc = ("BC_Topo_");
	std::string midname;
	std::string tailname = (".in");

	ss << ProcId;
	ss >> midname;

	std::string filename_bc = forename_bc + midname + tailname;

	std::ifstream myFile;
	myFile.open(filename_bc);

	if (myFile)
	{
		int nline=0;
		std::cout<<"Reading BC input file..."<<std::endl;
		while(readBCTopo(myFile,nline)&&(6!=nline))
		{
			int type=(BCTopo.rbegin())->Type;
			int id=(BCTopo.rbegin())->id;
			std::cout<<"Face "<<nline<<" type "<<BCTypes[type+1];
			if(type!=0)	std::cout<<" "<<BCnames[id-1];
			std::cout<<std::endl;
			nline++;
		};
		
		std::cout<<"Inlet conditions are given by:"<<std::endl;	
		for(auto iter=BCvalue.begin();iter!=BCvalue.end();iter++)
		{
			for(auto iter_value=iter->begin();iter_value!=iter->end();iter_value++)
				std::cout<<*iter_value<<" ";
			std::cout<<std::endl;		
		}

		if(6!=BCTopo.size())
		{
			std::cout<<"Error: BC file:"<<filename_bc<<" is incomplete"<<std::endl;
			std::cout<<"       BC number is less than 6:"<<BCTopo.size()<<std::endl;
			exit(-1); 
		}
	}
	else
	{
		std::cout<<"file "<<filename_bc<<" does not exist!"<<std::endl;
		exit(-1);
	}

	myFile.close();
	std::cout<<"BC Topo has been read!"<<std::endl;
	//myMPI.setTopo(BCTopo);


}


void nuc3d::boundaryCondition::updateBuffer_xi(VectorBuffer &myBf,VectorField &myFlux)
{
	if ((-1)==BCTopo[0].Type)
	{
		(this->*mySetter[BCTopo[0].id])(myBf,myFlux,0);
	}

	if ((-1)==BCTopo[1].Type)
	{
		(this->*mySetter[BCTopo[1].id])(myBf,myFlux,1);
	}
}

void nuc3d::boundaryCondition::updateBuffer_eta(VectorBuffer &myBf,VectorField &myFlux)
{
	if ((-1)==BCTopo[2].Type)
	{
		(this->*mySetter[BCTopo[2].id])(myBf,myFlux,2);
	}

	if ((-1)==BCTopo[3].Type)
	{
		(this->*mySetter[BCTopo[3].id])(myBf,myFlux,3);
	}
}

void nuc3d::boundaryCondition::updateBuffer_zeta(VectorBuffer &myBf,VectorField &myFlux)
{
	if ((-1)==BCTopo[4].Type)
	{
		(this->*mySetter[BCTopo[4].id])(myBf,myFlux,4);
	}

	if ((-1)==BCTopo[5].Type)
	{
		(this->*mySetter[BCTopo[5].id])(myBf,myFlux,5);
	}
}

std::ifstream& nuc3d::boundaryCondition::readBCTopo(std::ifstream& ios,int iface)
{
	std::string line;
	std::vector<double> BCvalue0(0);
	if(std::getline(ios,line))
	{
		int type;
		int id;
		double value; 
		std::istringstream values(line);
		values>>type>>id;
		BCTopo.push_back(faceBC(type,id));
		if((-1)==type&&(1==id))
		{
			while(values>>value)
			{
				BCvalue0.push_back(value);
			}

			if(BCvalue0.size()!=5)
			{
				std::cout<<"BC value number is less than 5:"<<BCvalue0.size()<<std::endl;
				exit(-1);
			}
			else
			{
				BCvalue.push_back(BCvalue0);
			}
		}
	}

	return ios;
}


void nuc3d::boundaryCondition::setBC(PDEData3d &myPDE,
		physicsModel &myPhyMod,
		EulerData3D &myFluxes)
{
	for(auto iter=BCTopo.begin();iter!=BCTopo.end();iter++)
	{
		int iface=static_cast<int>(iter-BCTopo.begin());
		if (((-1)==iter->Type)&&(1==iter->id))
		{
			setBC_Inlet(myPDE,myFluxes,myPhyMod,iface);
		}
	}
}

void nuc3d::boundaryCondition::setBC_Inlet(PDEData3d &myPDE,
		EulerData3D &myFluxes,
		physicsModel &myMod,
		int iface)
{
	double T;
	double e;
	double alpha;

	double rho=*(BCvalue[iface].begin());
	double u=*(BCvalue[iface].begin()+1);
	double v=*(BCvalue[iface].begin()+2);
	double w=*(BCvalue[iface].begin()+3);
	double p=*(BCvalue[iface].begin()+4);

	myMod.prim2conPoint(rho, u, v, w, p, T, e, alpha);

	int nx=myFluxes.nx;
	int ny=myFluxes.ny;
	int nz=myFluxes.nz;

	VectorField &W=myFluxes.getPrimatives();
	VectorField &W0=myFluxes.getAcoustics();

	int ibeg=0,iend=nx;
	int jbeg=0,jend=ny;
	int kbeg=0,kend=nz;

	switch (iface) {
		case 0:
			iend=1;
			break;
		case 1:
			ibeg=nx-1;
			break;
		case 2:
			jend=1;
			break;
		case 3:
			jbeg=ny-1;
			break;
		case 4:
			kend=1;
			break;
		case 5:
			kbeg=nz-1;
			break;
	}

	for (int k=kbeg; k<kend; k++)
	{
		for (int j=jbeg; j<jend; j++)
		{
			for (int i=ibeg; i<iend; i++)
			{
				W[0].setValue(i, j, k, rho);
				W[1].setValue(i, j, k, u);
				W[2].setValue(i, j, k, v);
				W[3].setValue(i, j, k, w);
				W[4].setValue(i, j, k, p);

				W0[0].setValue(i, j, k, T);
				W0[1].setValue(i, j, k, e);
				W0[2].setValue(i, j, k, alpha);
			}
		}
	}

}


void nuc3d::boundaryCondition::setBuffer_inlet(VectorBuffer &myBf,
		VectorField &myFlux,
		int iface)
{
	switch(iface)
	{
		case 0:
		case 1:
			setBuffer_inlet_xi(myBf, myFlux, iface);
			break;
		case 2:
		case 3:
			setBuffer_inlet_eta(myBf, myFlux, iface);
			break;
		case 4:
		case 5:
			setBuffer_inlet_zeta(myBf, myFlux, iface);
			break;
	}
}

void nuc3d::boundaryCondition::setBuffer_inlet_xi(VectorBuffer &myBf,
		VectorField &myFlux,
		int iface)
{
	auto beg=myFlux.begin();
	auto end=myFlux.end();

	int sizeX=myBf[0].BufferRecv[iface].getSizeX();
	int sizeY=myBf[0].BufferRecv[iface].getSizeY();
	int sizeZ=myBf[0].BufferRecv[iface].getSizeZ();

	int ix=(0==iface)?0:(beg->getSizeX()-1);

	for (auto iter=beg;iter!=end;iter++)
	{
		for (int k=0; k<sizeZ; k++)
		{
			for (int j=0; j<sizeY; j++)
			{
				for (int i=0; i<sizeX; i++)
				{
					double value=iter->getValue(ix, j, k);
					myBf[iter-beg].BufferRecv[iface].setValue(i, j, k, value);
				}
			}
		}
	}
}

void nuc3d::boundaryCondition::setBuffer_inlet_eta(VectorBuffer &myBf,
		VectorField &myFlux,
		int iface)
{
	auto beg=myFlux.begin();
	auto end=myFlux.end();

	int sizeX=myBf[0].BufferRecv[iface].getSizeX();
	int sizeY=myBf[0].BufferRecv[iface].getSizeY();
	int sizeZ=myBf[0].BufferRecv[iface].getSizeZ();

	int jx=(2==iface)?0:(beg->getSizeY()-1);

	for (auto iter=beg;iter!=end;iter++)
	{
		for (int k=0; k<sizeZ; k++)
		{
			for (int j=0; j<sizeY; j++)
			{
				for (int i=0; i<sizeX; i++)
				{
					double value=iter->getValue(i, jx, k);
					myBf[iter-beg].BufferRecv[iface].setValue(i, j, k, value);
				}
			}
		}
	}


}

void nuc3d::boundaryCondition::setBuffer_inlet_zeta(VectorBuffer &myBf,
		VectorField &myFlux,
		int iface)
{
	auto beg=myFlux.begin();
	auto end=myFlux.end();

	int sizeX=myBf[0].BufferRecv[iface].getSizeX();
	int sizeY=myBf[0].BufferRecv[iface].getSizeY();
	int sizeZ=myBf[0].BufferRecv[iface].getSizeZ();

	int kx=(4==iface)?0:(beg->getSizeZ()-1);

	for (auto iter=beg;iter!=end;iter++)
	{
		for (int k=0; k<sizeZ; k++)
		{
			for (int j=0; j<sizeY; j++)
			{
				for (int i=0; i<sizeX; i++)
				{
					double value=iter->getValue(i, j, kx);
					myBf[iter-beg].BufferRecv[iface].setValue(i, j, k, value);
				}
			}
		}
	}


}

void nuc3d::boundaryCondition::setBC_Outlet(PDEData3d &myPDE,
		EulerData3D &myFluxes,
		physicsModel &myMod,
		int iface)
{
	/* this interface is left for future usage */
}

void nuc3d::boundaryCondition::setBuffer_outlet(VectorBuffer &myBf,
		VectorField &myFlux,
		int iface)
{
	switch(iface)
	{
		case 0:
		case 1:
			setBuffer_outlet_xi(myBf, myFlux, iface);
			break;
		case 2:
		case 3:
			setBuffer_inlet_eta(myBf, myFlux, iface);
			break;
		case 4:
		case 5:
			setBuffer_inlet_zeta(myBf, myFlux, iface);
			break;
	}
}


void nuc3d::boundaryCondition::setBuffer_outlet_xi(VectorBuffer &myBf,
		VectorField &myFlux,
		int iface)
{
	auto beg=myFlux.begin();
	auto end=myFlux.end();

	int sizeX=myBf[0].BufferRecv[iface].getSizeX();
	int sizeY=myBf[0].BufferRecv[iface].getSizeY();
	int sizeZ=myBf[0].BufferRecv[iface].getSizeZ();

	int ibeg=(0==iface)?0:(beg->getSizeX()-1);
	int dir=(0==iface)?1:(-1);

	int ibeg_bf=(0==iface)?(sizeX-1):0;
	int iend_bf=(0==iface)?0:(sizeX-1);
	int dir_bf=(0==iface)?-1:1;

	for (auto iter=beg;iter!=end;iter++)
	{
		for (int k=0; k<sizeZ; k++)
		{
			for (int j=0; j<sizeY; j++)
			{
				double qn=iter->getValue(ibeg ,j, k);
				double qn_1=iter->getValue(ibeg+dir ,j, k);
				double k0=qn-qn_1/3.0;
				double q=qn*4.0/3.0-qn_1/3.0;

				for (int i=ibeg_bf; i!=(iend_bf+dir_bf); i+=dir_bf)
				{
					myBf[iter-beg].BufferRecv[iface].setValue(i, j, k, q);
					q=q/3.0+k0;
				}
			}
		}
	}


}

void nuc3d::boundaryCondition::setBuffer_outlet_eta(VectorBuffer &myBf,
		VectorField &myFlux,
		int iface)
{
	auto beg=myFlux.begin();
	auto end=myFlux.end();

	int sizeX=myBf[0].BufferRecv[iface].getSizeX();
	int sizeY=myBf[0].BufferRecv[iface].getSizeY();
	int sizeZ=myBf[0].BufferRecv[iface].getSizeZ();

	int jbeg=(2==iface)?0:(beg->getSizeY()-1);
	int dir=(2==iface)?1:(-1);

	int jbeg_bf=(2==iface)?(sizeY-1):0;
	int jend_bf=(2==iface)?0:(sizeY-1);
	int dir_bf=(2==iface)?-1:1;

	std::vector<double> q_vec(sizeX,0.0);
	std::vector<double> k_vec(sizeX,0.0);

	for (auto iter=beg;iter!=end;iter++)
	{
		for (int k=0; k!=sizeZ; k++)
		{
			for(auto iterQ=q_vec.begin();iterQ!=q_vec.end();iterQ++)
			{
				double qn=iter->getValue(static_cast<int>(iterQ-q_vec.begin()) ,jbeg, k);
				double qn_1=iter->getValue(static_cast<int>(iterQ-q_vec.begin()) ,jbeg+dir, k);
				double k0=qn-qn_1/3.0;
				double q=qn*4.0/3.0-qn_1/3.0;

				*iterQ=q;
				k_vec[iterQ-q_vec.begin()]=k0;
			}

			for (int j=jbeg_bf; j!=(jend_bf+dir_bf); j+=dir_bf)
			{
				for (int i=0; i!=sizeX; i++)
				{
					double q=q_vec[i];
					myBf[iter-beg].BufferRecv[iface].setValue(i, j, k, q);
					q_vec[i]=q/3.0+k_vec[i];
				}
			}
		}
	}

}

void nuc3d::boundaryCondition::setBuffer_outlet_zeta(VectorBuffer &myBf,
		VectorField &myFlux,
		int iface)
{
	auto beg=myFlux.begin();
	auto end=myFlux.end();

	int sizeX=myBf[0].BufferRecv[iface].getSizeX();
	int sizeY=myBf[0].BufferRecv[iface].getSizeY();
	int sizeZ=myBf[0].BufferRecv[iface].getSizeZ();

	int kbeg=(4==iface)?0:(beg->getSizeZ()-1);
	int dir=(4==iface)?1:(-1);

	int kbeg_bf=(4==iface)?(sizeZ-1):0;
	int kend_bf=(4==iface)?0:(sizeZ-1);
	int dir_bf=(4==iface)?-1:1;

	std::vector<double> q_vec(sizeX,0.0);
	std::vector<double> k_vec(sizeX,0.0);

	for (auto iter=beg;iter!=end;iter++)
	{
		for (int j=0; j!=sizeY; j++)
		{
			for(auto iterQ=q_vec.begin();iterQ!=q_vec.end();iterQ++)
			{
				double qn=iter->getValue(static_cast<int>(iterQ-q_vec.begin()) ,j, kbeg);
				double qn_1=iter->getValue(static_cast<int>(iterQ-q_vec.begin()) ,j, kbeg+dir);
				double k0=qn-qn_1/3.0;
				double q=qn*4.0/3.0-qn_1/3.0;

				*iterQ=q;
				k_vec[iterQ-q_vec.begin()]=k0;
			}

			for (int k=kbeg_bf; k!=(kend_bf+dir_bf); k+=dir_bf)
			{
				for (int i=0; i!=sizeX; i++)
				{
					double q=q_vec[i];
					myBf[iter-beg].BufferRecv[iface].setValue(i, j, k, q);
					q_vec[i]=q/3.0+k_vec[i];
				}
			}
		}
	}


}

void nuc3d::boundaryCondition::setBC_wall(PDEData3d &,
		EulerData3D &,
		physicsModel &myMod,
		int)
{
	/* this interface is left for future usage */
}

void nuc3d::boundaryCondition::setBuffer_wall(VectorBuffer &myBf,
		VectorField &myFlux,
		int iface)
{
	switch(iface)
	{
		case 0:
		case 1:
			setBuffer_wall_xi(myBf, myFlux, iface);
			break;
		case 2:
		case 3:
			setBuffer_wall_eta(myBf, myFlux, iface);
			break;
		case 4:
		case 5:
			setBuffer_wall_zeta(myBf, myFlux, iface);
			break;
	}

}

void nuc3d::boundaryCondition::setBuffer_wall_xi(VectorBuffer &myBf,
		VectorField &myFlux,
		int iface)
{
	std::vector<double> q_vec(3,0.0);
	auto beg=myFlux.begin();
	auto end=myFlux.end();

	int sizeX=myBf[0].BufferRecv[iface].getSizeX();
	int sizeY=myBf[0].BufferRecv[iface].getSizeY();
	int sizeZ=myBf[0].BufferRecv[iface].getSizeZ();

	int ibeg=(0==iface)?0:(beg->getSizeX()-1);
	int dir=(0==iface)?1:(-1);

	int ibeg_bf=(0==iface)?(sizeX-1):0;
	int iend_bf=(0==iface)?0:(sizeX-1);
	int dir_bf=(0==iface)?-1:1;

	for (auto iter=beg;iter!=end;iter++)
	{
		for (int k=0; k<sizeZ; k++)
		{
			for (int j=0; j<sizeY; j++)
			{
				for(int i=0;i!=sizeX;i++)
				{
					q_vec[i]=iter->getValue(ibeg+i*dir,j, k);
				}

				for (int i=ibeg_bf; i!=(iend_bf+dir_bf); i+=dir_bf)
				{
					myBf[iter-beg].BufferRecv[iface].setValue(i, j, k, q_vec[i-ibeg_bf]);
				}
			}
		}
	}
}

void nuc3d::boundaryCondition::setBuffer_wall_eta(VectorBuffer &myBf,
		VectorField &myFlux,
		int iface)
{
	auto beg=myFlux.begin();
	auto end=myFlux.end();

	int sizeX=myBf[0].BufferRecv[iface].getSizeX();
	int sizeY=myBf[0].BufferRecv[iface].getSizeY();
	int sizeZ=myBf[0].BufferRecv[iface].getSizeZ();

	int jbeg=(2==iface)?0:(beg->getSizeY()-1);

	int jbeg_bf=(2==iface)?(sizeY-1):0;
	int jend_bf=(2==iface)?0:(sizeY-1);
	int dir_bf=(2==iface)?-1:1;

	for (auto iter=beg;iter!=end;iter++)
	{
		for (int k=0; k<sizeZ; k++)
		{
			for (int j=jbeg_bf; j!=(jend_bf+dir_bf); j+=dir_bf)
			{
				for (int i=0; i!=sizeX; i++)
				{
					double q=iter->getValue(i, jbeg-(j-jbeg), k);
					myBf[iter-beg].BufferRecv[iface].setValue(i, j, k, q);
				}
			}
		}
	}

}

void nuc3d::boundaryCondition::setBuffer_wall_zeta(VectorBuffer &myBf,
		VectorField &myFlux,
		int iface)
{
	auto beg=myFlux.begin();
	auto end=myFlux.end();

	int sizeX=myBf[0].BufferRecv[iface].getSizeX();
	int sizeY=myBf[0].BufferRecv[iface].getSizeY();
	int sizeZ=myBf[0].BufferRecv[iface].getSizeZ();

	int kbeg=(4==iface)?0:(beg->getSizeZ()-1);

	int kbeg_bf=(4==iface)?(sizeZ-1):0;
	int kend_bf=(4==iface)?0:(sizeZ-1);
	int dir_bf=(4==iface)?-1:1;

	for (auto iter=beg;iter!=end;iter++)
	{
		for (int k=kbeg_bf; k!=(kend_bf+dir_bf); k+=dir_bf)
		{
			for (int j=0; j!=sizeY; j++)
			{
				for (int i=0; i!=sizeX; i++)
				{
					double q=iter->getValue(i, j, kbeg-(k-kbeg_bf));
					myBf[iter-beg].BufferRecv[iface].setValue(i, j, k, q);
				}
			}
		}
	}

}


void nuc3d::boundaryCondition::setBuffer_symm(VectorBuffer &,VectorField &,int)
{

}
void nuc3d::boundaryCondition::setBuffer_period(VectorBuffer &,VectorField &,int)
{

}
