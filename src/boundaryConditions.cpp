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
        while(readBCTopo(myFile));
        if(6!=BCTopo.size())
        {
            std::cout<<"Error: BC file:"<<filename_bc<<" is incomplete"<<std::endl;
            std::cout<<"       BC number is less than 6"<<std::endl;
            
        }
        
        for(auto iter=BCTopo.begin();iter!=BCTopo.end();iter++)
        {
            int iface=static_cast<int>(iter-BCTopo.begin());
            if(1==iter->id)
            {
                readBCValue(myFile,iface);
            }
        }
    }
    else
    {
        std::cout<<"file "<<filename_bc<<" does not exist!"<<std::endl;
        exit(-1);
    }
    
    myFile.close();
    myMPI.setTopo(BCTopo);
    
}


void nuc3d::boundaryCondition::updateBC(VectorBuffer &,EulerData3D &)
{
    
}

std::istream& nuc3d::boundaryCondition::readBCTopo(std::istream& ios)
{
    double value0;
    double value1;
    
    ios >> value0 >> value1;
    
    BCTopo.push_back(faceBC(value0,value1));
    
    return ios;
}

std::istream& nuc3d::boundaryCondition::readBCValue(std::istream &ios,int iface)
{
    std::string line;
    double value;
    std::vector<double> BCvalue0;
    if(std::getline(ios,line))
    {
        std::istringstream values(line);
        while(values>>value)
        {
            BCvalue0.push_back(value);
        }
        
        if(BCvalue0.size()<5)
        {
            std::cout<<"BC value number is less than 5!"<<std::endl;
            exit(-1);
            
        }
    }
    else
    {
        std::cout<<"BC file is incomplete!"<<std::endl;
        exit(-1);
    }
    
    BCvalue[iface]=BCvalue0;
    return ios;
    
}

void nuc3d::boundaryCondition::setBC(PDEData3d &myPDE,
                                     physicsModel &myPhyMod,
                                     EulerData3D &myFluxes)
{
    for(auto iter=BCTopo.begin();iter!=BCTopo.end();iter++)
    {
        int iface=static_cast<int>(iter-BCTopo.begin());
        if (((-1)==iter->Type)&&(1==iter->id)) {
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
    
    
    
}
