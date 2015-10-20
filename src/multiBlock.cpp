#ifndef meshBlock_cpp
#define meshBlock_cpp
#include "multiBlock.h"
/**************************************************************************************
 Definition of class communicator: base class for communicators
 **************************************************************************************/
/**************************************************************************************
 Definition of constructors and destructors
 **************************************************************************************/
nuc3d::MeshBlock::MeshBlock(int& argc, char **& argv):
myComm(argc, argv),
myCtrler(),
myPhys()
{
    
}

nuc3d::MeshBlock::MeshBlock(int nx0,int ny0,int nz0,
                            int bfwidth,
                            int blkId,
                            const physicsModel &PhysMod,
                            const fieldOperator3d &FieldOp,
                            const MPIComunicator3d_nonblocking &Commtor):
nx(nx0),ny(ny0),nz(nz0),
bufferWidth(FieldOp.bufferSize),
myBlkId(blkId),
xyz(3,Field(nx0+1,ny0+1,nz0+1)),
xyz_center(3,Field(nx0,ny0,nz0))
{
    
    myPhys=PhysMod;
    myOperator=FieldOp;
    myComm=Commtor;
    
    if("Euler3d"==myPhysBLK.getMyModelName())
    {
        myEuler=std::make_shared<EulerData3D>(nx0,ny0,nz0,neqs);
    }
    else if ("EulerReactive3d"==myPhysBLK.getMyModelName())
    {
        myEuler=std::make_shared<EulerReactiveData3D>(nx0,ny0,nz0,neqs);
    }
    else if ("NaiverStokes3d"==myPhysBLK.getMyModelName())
    {
        myEuler=std::make_shared<NaiverStokesData3d>(nx0,ny0,nz0,neqs);
    }
    else if ("NaiverStokesReactive3d"==myPhysBLK.getMyModelName())
    {
        myEuler=std::make_shared<NaiverStokesReactiveData3d>(nx0,ny0,nz0,neqs);
    }
    else
    {
        std::cout <<"Model Name:"<<myPhysBLK.getMyModelName()
        <<"does not exist!"
        <<std::endl;
        exit(-1);
    }
    
    for(int i=0;i<myPhysBLK.getEqNum();i++)
        myBuffers.push_back(bufferData(nx0,ny0,nz0,bfwidth));
};

nuc3d::MeshBlock::~MeshBlock()
{
    
};
/**************************************************************************************
 Definition of member functions
 **************************************************************************************/
void nuc3d::MeshBlock::input()
{
    std::stringstream ss;
    int ProcId = myComm.getMyId();
    int nx0, ny0, nz0;
    std::string forename_mesh = (".mesh_");
    std::string midname;
    std::string tailname = (".dat");
    
    ss << ProcId;
    ss >> midname;
    
    std::string filename_mesh = forename + midname + tailname;
    
    
    std::ifstream myFile(filename);
    if (myFile)
    {
        myFile >> nx >> ny >> nz;
        myPDE＝PDEData3d(nx,ny,nz,myPhys.getEqNum());
        
        myOperator=fieldOperator3d(myPDE.getQ());
        myBlkId=myComm.getMyId();
        
        bufferWidth=myOperator.getBufferSize();
        
        for(int i=0;i<3;i++)
        {
            xyz.push_back(Field(nx+1,ny+1,nz+1));
        }
        
        for(int i=0;i<3;i++)
            xyz_center.push_back(Field(nx,ny,nz));
        
        if("Euler3d"==myPhysBLK.getMyModelName())
        {
            myEuler=std::make_shared<EulerData3D>(nx,ny,nz,myPhys.getEqNum());
        }
        else if ("EulerReactive3d"==myPhysBLK.getMyModelName())
        {
            myEuler=std::make_shared<EulerReactiveData3D>(nx,ny,nz,myPhys.getEqNum());
        }
        else if ("NaiverStokes3d"==myPhysBLK.getMyModelName())
        {
            myEuler=std::make_shared<NaiverStokesData3d>(nx,ny,nz,myPhys.getEqNum());
        }
        else if ("NaiverStokesReactive3d"==myPhysBLK.getMyModelName())
        {
            myEuler=std::make_shared<NaiverStokesReactiveData3d>(nx,ny,nz,myPhys.getEqNum());
        }
        else
        {
            std::cout <<"Model Name:"<<myPhysBLK.getMyModelName()
            <<"does not exist!"
            <<std::endl;
            exit(-1);
        }
        
        for(int i=0;i<myPhysBLK.getEqNum();i++)
            myBuffers.push_back(bufferData(nx0,ny0,nz0,bfwidth));
        
        double x,y,z;
        
        for(int k=0;k<nz+1;k++)
        {
            for(int j=0;j<ny+1;j++)
            {
                for(int i=0;i<nx+1;i++)
                {
                    myFile>>x>>y>>z;
                    xyz[0].setValue(i, j, k, x);
                    xyz[1].setValue(i, j, k, y);
                    xyz[2].setValue(i, j, k, z);
                }
            }
        }
    }
    else
    {
        std::cout << "file " << filename << " does not exist" << std::endl;
        myComm.AbortMPI();
    }
    myFile.close();
    readData(myFile,xyz);
    readData(myFile,myFluxes->W_Euler);
}

void nuc3d::MeshBlock::inputXYZ()
{
    std::stringstream ss;
    int ProcId = myComm.getMyId();
    int nx0, ny0, nz0;
    std::string forename = ("mesh_");
    std::string midname;
    std::string tailname = (".dat");
    
    ss << ProcId;
    ss >> midname;
    
    std::string filename = forename + midname + tailname;
    
    
    std::ifstream myFile(filename);
    if (myFile)
    {
        myFile >> nx0 >> ny0 >> nz0
        double x,y,z;
        for(int k=0;k<nz0;k++)
            for(int j=0;j<ny0;j++)
                for(int i=0;i<nx0;i++)
                {
                    myFile>>x>>y>>z;
                    xyz[0].setValue(i, j, k, x);
                    xyz[1].setValue(i, j, k, y);
                    xyz[2].setValue(i, j, k, z);
                }
    }
    else
    {
        std::cout << "file " << filename << " does not exist" << std::endl;
        myComm.AbortMPI();
    }
    
    myFile.close();
    
}

void nuc3d::MeshBlock::readData(std::ifstream &myFile, VectorField &dataVec)
{
    double value;
    
    int nx=(dataVec.begin())->getSizeX();
    int ny=(dataVec.begin())->getSizeY();
    int nz=(dataVec.begin())->getSizeZ();
    
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
    int nx=(dataVec.begin())->getSizeX();
    int ny=(dataVec.begin())->getSizeY();
    int nz=(dataVec.begin())->getSizeZ();
    
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
    {
        //finish this part at oct 19 2015
    }
    
}


void nuc3d::MeshBlock::initialPDE()
{
    myPhysBLK.initial(myPDE, myFluxes);
}

void nuc3d::MeshBlock::solve()
{
    solveRiemann();
    solveInvicidFlux();
    solveViscousFLux();
    solveGetRHS();
    solveIntegral();
    
}

void nuc3d::MeshBlock::solveRiemann()
{
    myPhysBLK.solve(myPDE, myFluxes);
}

void nuc3d::MeshBlock::solveInvicidFlux()
{
}

void nuc3d::MeshBlock::solveViscousFLux()
{
    myFluxes->solveLocal();
}

void nuc3d::MeshBlock::solveGetRHS()
{
    
}

void nuc3d::MeshBlock::solveIntegral()
{
}


void nuc3d::MeshBlock::output(std::ofstream &)
{
    
}
/**************************************************************************************
 End of definition
 **************************************************************************************/
#endif