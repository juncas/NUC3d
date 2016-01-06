//
//  singleBlock.hpp
//  NUC3d
//
//  Created by Jun Peng on 15/10/20.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//
#ifndef singleBlock_hpp
#define singleBlock_hpp
#include <cstdlib>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include "field.h"
#include "MPICommunicator.h"
#include "block.h"
#include "physicsModel.h"
#include "IOcontroller.h"
#include "fieldOperator.h"

namespace nuc3d
{
    class MPIComunicator3d_nonblocking;
    class IOController;
    class physicsModel;
    class fieldOperator3d;
    class boundaryCondition;
    class block;
    
    class singleBlock
    {
        MPIComunicator3d_nonblocking myMPI;
        IOController myIO;
        physicsModel myPhys;
        fieldOperator3d myOperator;
        boundaryCondition myBC;
        
        block myBlock;
                
    public:
        singleBlock();
        
        ~singleBlock();
        
        void loop();
        
    private:
        
        void initialBlock();
        
        void solvePDE();
        
        void postprocess();
        
    private:
        void readData(std::ifstream &, VectorField &);
        void writeData(std::ofstream &, VectorField &);
                        
    };
}
#endif /* singleBlock_hpp */
