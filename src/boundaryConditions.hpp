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
#include <stdio.h>
#include "field.h"
#include "mpi.h"
#include "bufferData.hpp"

namespace nuc3d
{
    class EulerData3D;
    
    class BufferMetric
    {
        Field jacobian;
        
        VectorField xi_xyz;
        VectorField eta_xyz;
        VectorField zeta_xyz;
    
    public:
        BufferMetric(int,int,int);
        ~BufferMetric();
    };
    

    class boundaryCondition
    {
        int myBC[6];
        int myNeibBlk[6];
        int myNeiBlkFace[6];
        
        std::vector<double> BCvalue;
        
        std::vector<std::shared_ptr<BufferMetric>> myBufferGeo;
    public:
        boundaryCondition();
        ~boundaryCondition();
    
    public:
        void setBC(VectorBuffer &,EulerData3D &);
        void applyBC(VectorBuffer &,EulerData3D &);
        void updateBC(VectorBuffer &,EulerData3D &);
        void initialBC(VectorBuffer &);
    };
}
#endif /* boundaryConditions_hpp */
