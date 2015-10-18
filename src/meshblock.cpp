#ifndef meshBlock_cpp
#define meshBlock_cpp
#include "meshblock.h"
/**************************************************************************************
 Definition of class communicator: base class for communicators
 **************************************************************************************/
/**************************************************************************************
 Definition of constructors and destructors
 **************************************************************************************/
nuc3d::MeshBlock::MeshBlock(int nx0,int ny0,int nz0,int bfwidth,const physicsModel &myPhys):
nx(nx0),ny(ny0),nz(nz0),
bufferWidth(bfwidth),
xyz(3,Field(nx0,ny0,nz0)),
myPDE(nx0,ny0,nz0,myPhys.getEqNum()),
myPhysBLK(myPhys)
{
    myOperator(myPDE.Q_Euler);
    
    if("Euler"==myPhysBLK.getMyModelName())
    {
        myEuler=std::make_shared<EulerData3D>(nx0,ny0,nz0,neqs);
    }
    else if ("EulerReactive"==myPhysBLK.getMyModelName())
    {
        myEuler=std::make_shared<EulerReactiveData3D>(nx0,ny0,nz0,neqs);
    }
    else if ("NaiverStokes"==myPhysBLK.getMyModelName())
    {
        myEuler=std::make_shared<NaiverStokesData3d>(nx0,ny0,nz0,neqs);
    }
    else if ("NaiverStokesReactive"==myPhysBLK.getMyModelName())
    {
        myEuler=std::make_shared<NaiverStokesReactiveData3d>(nx0,ny0,nz0,neqs);
    }
    else
    {
        std::cout <<"Model Name:"<<myPhysBLK.getMyModelName()
        <<"does not exist!"
        <<std::endl;
    }
    
    for(int i=0;i<neqs;i++)
        myBuffers.push_back(bufferData(nx0,ny0,nz0,bfwidth));
};

nuc3d::MeshBlock::~MeshBlock()
{
    
};
/**************************************************************************************
 Definition of member functions
 **************************************************************************************/
void nuc3d::MeshBlock::input(std::ifstream &myFile)
{
    readData(myFile,xyz);
    readData(myFile,myEuler.W_Euler);
}

void nuc3d::MeshBlock::readData(std::ifstream &myFile, VectorField &dataVec)
{
    double value;
    
    for (auto iter = dataVec.begin(); iter != dataVec.end(); ++iter)
    {
        for(int k=0;k<nz;k++)
            for(int j=0;j<ny;j++)
                for(int i=0;i<nx;i++)
                {
                    myFile>>value;
                    iter->setValue(i, j, k, value);
                }
    }
    
}

void nuc3d::MeshBlock::writeData(std::ofstream &myFile, VectorField &dataVec)
{
    double value;
    
    for (auto iter = dataVec.begin(); iter != dataVec.end(); ++iter)
    {
        for(int k=0;k<nz;k++)
            for(int j=0;j<ny;j++)
                for(int i=0;i<nx;i++)
                {
                    value=iter->getValue(i,j,k);
                    myFile<<value<<std::endl;
                }
    }
    
}

void nuc3d::MeshBlock::initial()
{
    initialXYZ();
}

void nuc3d::MeshBlock::initialXYZ()
{
    for(auto iterXYZ=xyz.begin();iterXYZ!=xyz.end();++iterXYZ)
        
        }

void nuc3d::MeshBlock::solve()
{
    
}

void nuc3d::MeshBlock::output(std::ofstream &)
{
    
}
/**************************************************************************************
 End of definition
 **************************************************************************************/
#endif