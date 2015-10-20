//
//  singleBlock.hpp
//  NUC3d
//
//  Created by Jun Peng on 15/10/20.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#ifndef singleBlock_hpp
#define singleBlock_hpp
#include "block.hpp"

namespace nuc3d
{
    
    class singleBlock : public block
    {
    public:
        singleBlock();
        ~singleBlock();
        
        void solveRiemann();
        void solveBoundaryConditions();
        void solveInvicidFlux();
        void solveViscousFLux();
        void solveGetRHS();
        void solveIntegral();
        
        void input(std::ifstream &);
        void initial();
        
        void output(std::ofstream &);
        
        void setMyId(int Id) { myBlkId = Id; };
        
    private:
        void readData(std::ifstream &, VectorField &);
        void writeData(std::ofstream &, VectorField &);
        
        void initialXYZ();
        void initialPDE();
        
        void putXYZ();
        void putEuler();
        
    };
}
#endif /* singleBlock_hpp */
