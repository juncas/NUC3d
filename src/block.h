//
//  block.hpp
//  NUC3d
//
//  Created by Jun Peng on 15/10/20.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//
#ifndef block_hpp
#define block_hpp
#include <cstdlib>
#include <memory>
#include "field.h"
#include "PDEData3d.hpp"
#include "bufferData.hpp"

namespace nuc3d
{
    class PDEData3d;
    class EulerData3D;
    class bufferData;
    class physicsModel;
    class MPIComunicator3d_nonblocking;
    class boundaryCondition;
    class IOController;
    
    class block
    {
    protected:
        int nx;
        int ny;
        int nz;
	int bfsize;                
        VectorField xyz;
        VectorField xyz_center;
        
        PDEData3d myPDE;
        
        /*
         this shared point could be:
         - EulerData3D;
         - EulerReactiveData3D;
         - NaiverStokesData3d;
         - NaiverStokesReactiveData3d;
         */
        std::shared_ptr<EulerData3D> myFluxes;
        
        VectorBuffer mybuffer;
        
        double time;
        double dt;
        int istep;
        double RES;

    public:
        block();
        ~block();
        void initial(fieldOperator3d &,
                     physicsModel &,
                     MPIComunicator3d_nonblocking &,
                     boundaryCondition &,
                     IOController &);
        void solve(fieldOperator3d &,
                   physicsModel &,
                   MPIComunicator3d_nonblocking &,
                   boundaryCondition &,
                   IOController &);
        void printStatus();
    private:
        void initialData(int,int,int,physicsModel &);
        void ReadXYZ();
        void ReadPDE();

    };
    
}
#endif /* block_hpp */
