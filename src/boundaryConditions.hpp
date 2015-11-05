//
//  boundaryConditions.hpp
//  NUC3d
//
//  Created by Jun Peng on 15/11/5.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#ifndef boundaryConditions_hpp
#define boundaryConditions_hpp
#include <memory>
#include <map>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include "field.h"
#include "mpi.h"
#include "bufferData.hpp"

namespace nuc3d
{
    class EulerData3D;
    class MPIComunicator3d_nonblocking;
    class PDEData3d;
    class physicsModel;
        
    class faceBC
    {
    public:
        int Type; // -1:= bc  n:= neighbour block id
        int id; // bc id or neighbour face id
        
        faceBC(int myType,int myID):Type(myType),id(myID){};
        ~faceBC();
    };

    class boundaryCondition
    {
        std::vector<faceBC> BCTopo;
        
        std::map<int,std::vector<double>> BCvalue;//(rho,rhou,rhow,rhoe)
        
    public:
        boundaryCondition();
        ~boundaryCondition();
    
    public:
        void setBC(PDEData3d &,
                   physicsModel &,
                   EulerData3D &);
        
        void updateBC(VectorBuffer &,EulerData3D &);
        
        void initialBC(VectorBuffer &,
                       MPIComunicator3d_nonblocking &);
    private:
        std::istream& readBCTopo(std::istream&);
        std::istream& readBCValue(std::istream&,int);
        
        void setBC_Inlet(PDEData3d &,
                         EulerData3D &,
                         physicsModel &myMod,
                         int);
        void setBC_Outlet();
        void setBC_Wall();
        void setBC_Periodic();
        void setBC_Symmetric();
    };
}

#endif /* boundaryConditions_hpp */
