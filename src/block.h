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
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <memory>
#include "field.h"
#include "physicsModel.h"

namespace nuc3d
{
    
    class block
    {
        friend class physicsModel;
    protected:
        int nx;
        int ny;
        int nz;
                
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
        
        block();
        void initial(int,int,int,physicsModel &);
        ~block();
    };
    
}
#endif /* block_hpp */
