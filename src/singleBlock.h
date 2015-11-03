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
        
        void initialPDE();
        
        void solvePDE();
        
        void output();
        
	private:

		void solveInvicidFluxL(EulerFlux &, std::vector<bufferData> &, int);

		void solveInvicidFluxR(EulerFlux &, std::vector<bufferData> &, int);

        void solveRiemann();
        
        void solveBoundaryConditions();
        
        void solveInvicidFlux();
        
        void solveViscousFLux(EulerData3D &myEuler);
        void solveViscousFLux(EulerReactiveData3D &myEuler);
        void solveViscousFLux(NaiverStokesData3d &myEuler);
        void solveViscousFLux(NaiverStokesReactiveData3d &myEule);
        
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
