//
//  singleBlock.hpp
//  NUC3d
//
//  Created by Jun Peng on 15/10/20.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//
#ifndef singleBlock_hpp
#define singleBlock_hpp
#include "block.h"
#include "physicsModel.h"
#include "IOcontroller.h"
#include "fieldOperator.h"
#include "MPICommunicator.h"

namespace nuc3d
{
    
    class singleBlock : public block
    {
        MPIComunicator3d_nonblocking myComm;
        IOController myCtrler;
        physicsModel myPhys;
        fieldOperator3d myOperator;
    
    public:
        singleBlock(int&, char **&);
        ~singleBlock();
        
        void initialBlock();
        
        void initial();
        
        void solve();
        
        void output();
        
    private:
        void solveRiemann();
        
        void solveBoundaryConditions();
        
        void solveInvicidFlux();
        
        void solveViscousFLux();
        
        void solveGetRHS();
        
        void solveIntegral();
        
        void printRES();
    
    private:
        void readData(std::ifstream &, VectorField &);
        void writeData(std::ofstream &, VectorField &);
        
        void initialXYZ();
        void initialPDE();
                
    };
}
#endif /* singleBlock_hpp */
