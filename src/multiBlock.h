#ifndef meshBlock_h
#define meshBlock_h
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
       
    class multiBlock : public block //contains field datas of a mesh block
    {
        friend class physicsModel;
        friend class meshGrid;
    protected:
        int myBlkId;
        
        physicsModel &myPhys;
        fieldOperator3d &myOperator;
        MPIComunicator3d_nonblocking &myComm;
        IOController &myCtrler;
        
        
        
        
        std::vector<bufferData> myBuffers;
        
    public:
        // for single block configuration
        multiBlock(int& argc, char **& argv);
        
        // for single block configuration
        multiBlock(int,int,int,int,int,
                  const physicsModel &,
                  const fieldOperator3d &,
                  const MPIComunicator3d_nonblocking &);
        
        ~MeshBlock();
        
    public:
        
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

#endif