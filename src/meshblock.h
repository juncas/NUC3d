#ifndef meshBlock_h
#define meshBlock_h

#include <cstdlib>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include "euler3d.h"
#include "MPICommunicator.h"

namespace nuc3d
{
    class MeshBlock//contains field datas of a mesh block
    {
        friend class physicsModel;
        friend class meshGrid;
        
        int myBlkId;
        
        int nx;
        int ny;
        int nz;
        
        int bufferWidth;
        
        VectorField xyz;
        VectorField xyz_center;
        
        PDEData3d myPDE;
        physicsModel &myPhysBLK;
        
        /*
         this shared point could be:
         - EulerData3D;
         - EulerReactiveData3D;
         - NaiverStokesData3d;
         - NaiverStokesReactiveData3d;
         */
        
        std::shared_ptr<EulerData3D> myFluxes;
        
        std::vector<bufferData> myBuffers;//finished
        
        fieldOperator3d myOperator;
        
        //Boundary setBoundary;
        //Boundary conditions should be owned by each grid
        
    public:
        
        MeshBlock(int,int,int,int,int,const physicsModel &);
        ~MeshBlock();
        
    public:
        
        void input(std::ifstream &);
        void initial();
        
        void output(std::ofstream &);
        
        void setMyId(int Id) { myBlkId = Id; };
        
    private:
        void readData(std::ifstream &, VectorField &);
        void writeData(std::ofstream &, VectorField &);
        
        void initialXYZ();
        
        void putXYZ();
        void putEuler();
    };
}

#endif