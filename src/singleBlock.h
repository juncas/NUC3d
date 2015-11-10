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
#include "block.h"
#include "physicsModel.h"
#include "IOcontroller.h"
#include "fieldOperator.h"
#include "MPICommunicator.h"

namespace nuc3d
{
    class boundaryCondition;
    
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
        
        void initialBlock();
        
        void solvePDE();
        
        void postprocess();
        
        void output();
        
	private:
        void solveRiemann();
        
        void solveBoundaryConditions();
        
        void solveInvicidFlux();
        
        void solveInvicidFluxL(EulerFlux &, std::vector<bufferData> &, int);
        
        void solveInvicidFluxR(EulerFlux &, std::vector<bufferData> &, int);
        
        void solveViscousFLux();
        
        void solveGetRHS();
        
        void solveIntegral(int step);
        
        void printRES();
    
    private:
        void readData(std::ifstream &, VectorField &);
        void writeData(std::ofstream &, VectorField &);
        
        void initialXYZ();
        void initialPDE();
                
    };
}
#endif /* singleBlock_hpp */
