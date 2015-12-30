#ifndef euler3d_cpp
#define euler3d_cpp
#include <cmath>
#include "euler3d.h"
#include "physicsModel.h"
#include "PDEData3d.hpp"
#include "bufferData.hpp"

/**************************************************************************************
 Member functions of class: EulerFlux
 **************************************************************************************/
nuc3d::EulerFlux::EulerFlux(int nx0,int ny0,int nz0,int neqs):
FluxL(neqs,Field(nx0,ny0,nz0)),
FluxR(neqs,Field(nx0,ny0,nz0)),
reconstFluxL(neqs,Field(nx0+1,ny0,nz0)),
reconstFluxR(neqs,Field(nx0+1,ny0,nz0)),
reconstFlux(neqs,Field(nx0+1,ny0,nz0)),
maxEigen(1.0)
{
    
}

void nuc3d::EulerFlux::combineFluxLR()
{
    for(auto iter=reconstFlux.begin();iter!=reconstFlux.end();iter++)
        *iter=reconstFluxL[iter-reconstFlux.begin()]+reconstFluxR[iter-reconstFlux.begin()];
}

nuc3d::EulerFlux::~EulerFlux()
{}

/**************************************************************************************
 Member functions of class: EulerData3D
 **************************************************************************************/
nuc3d::EulerData3D::EulerData3D( int nx0, int ny0, int nz0, int neqs):
nx(nx0),
ny(ny0),
nz(nz0),
jacobian(Field(nx0,ny0,nz0)),
xi_xyz(3,Field(nx0,ny0,nz0)),
eta_xyz(3,Field(nx0,ny0,nz0)),
zeta_xyz(3,Field(nx0,ny0,nz0)),
dx(3,Field(nx0,ny0,nz0)),
dy(3,Field(nx0,ny0,nz0)),
dz(3,Field(nx0,ny0,nz0)),
W0_Euler(3,Field(nx0,ny0,nz0)),
W_Euler(neqs,Field(nx0,ny0,nz0)),
Flux_xi(nx0,ny0,nz0,neqs),
Flux_eta(ny0,nz0,nx0,neqs),
Flux_zeta(nz0,nx0,ny0,neqs),
dfdxi(neqs,Field(nx0,ny0,nz0)),
dgdeta(neqs,Field(ny0,nz0,nx0)),
dhdzeta(neqs,Field(nz0,nx0,ny0)),
dt(0.0)
{
}

void nuc3d::EulerData3D::initialxyz(VectorField &XYZ_c)
{
    
}

nuc3d::VectorField& nuc3d::EulerData3D::getDrivativeXi()
{
    return this->dfdxi;
}

nuc3d::VectorField& nuc3d::EulerData3D::getDrivativeEta()
{
    return this->dgdeta;
}

nuc3d::VectorField& nuc3d::EulerData3D::getDrivativeZeta()
{
    return this->dhdzeta;
}

nuc3d::EulerFlux& nuc3d::EulerData3D::getFluxXi()
{
    return Flux_xi;
}

nuc3d::EulerFlux& nuc3d::EulerData3D::getFluxEta()
{
    return Flux_eta;
}

nuc3d::EulerFlux& nuc3d::EulerData3D::getFluxZeta()
{
    return Flux_zeta;
}

nuc3d::VectorField& nuc3d::EulerData3D::getPrimatives()
{
    return W_Euler;
}


nuc3d::VectorField& nuc3d::EulerData3D::getAcoustics()
{
    return W0_Euler;
}


void nuc3d::EulerData3D::setDerivativesInv()
{
    Flux_xi.combineFluxLR();
    Flux_eta.combineFluxLR();
    Flux_zeta.combineFluxLR();
    
    EulerData3D::setDerivativesXiInv();
    EulerData3D::setDerivativesEtaInv();
    EulerData3D::setDerivativesZetaInv();
}

void nuc3d::EulerData3D::setDerivativesXiInv()
{
    VectorField &flux=Flux_xi.reconstFlux;
    for (auto iter=flux.begin() ; iter!=flux.end(); iter++)
    {
        double *f=iter->getDataPtr();
        double *df=dfdxi[iter-flux.begin()].getDataPtr();
        
        int nx0=iter->getSizeX();
        int ny0=iter->getSizeY();
        int nz0=iter->getSizeZ();
        
        int nx1=dfdxi[iter-flux.begin()].getSizeX();
        int ny1=dfdxi[iter-flux.begin()].getSizeY();
        int nz1=dfdxi[iter-flux.begin()].getSizeZ();
        
        for (int k=0; k<nz0; k++)
        {
            for (int j=0; j<ny0; j++)
            {
                for (int i=1; i<nx0; i++)
                {
                    int idx=nx0*ny0*k+nx0*j+i;
                    int idx1=nx1*ny1*k+nx1*j+i;
                    double df_loc=f[idx]-f[idx-1];
                    
                    df[idx1-1]=df_loc;
                }
            }
        }
    }
}

void nuc3d::EulerData3D::setDerivativesEtaInv()
{
    VectorField &flux=Flux_eta.reconstFlux;
    
    for (auto iter=flux.begin() ; iter!=flux.end(); iter++)
    {
        double *f=iter->getDataPtr();
        int nx0=iter->getSizeX();
        int ny0=iter->getSizeY();
        int nz0=iter->getSizeZ();

        double *dg=dgdeta[iter-flux.begin()].getDataPtr();
        int nx1=dgdeta[iter-flux.begin()].getSizeX();
        int ny1=dgdeta[iter-flux.begin()].getSizeY();
        int nz1=dgdeta[iter-flux.begin()].getSizeZ();
        
        for (int k=0; k<nz0; k++)
        {
            for (int j=0; j<ny0; j++)
            {
                for (int i=1; i<nx0; i++)
                {
                    int idx=nx0*ny0*k+nx0*j+i;
                    int idx1=nx1*ny1*k+nx1*j+i;
                    double df_loc=f[idx]-f[idx-1];
                    
                    dg[idx1-1]=df_loc;
                }
            }
        }
    }
}

void nuc3d::EulerData3D::setDerivativesZetaInv()
{
    VectorField &flux=Flux_zeta.reconstFlux;
    
    for (auto iter=flux.begin() ; iter!=flux.end(); iter++)
    {
        double *f=iter->getDataPtr();
        int nx0=iter->getSizeX();
        int ny0=iter->getSizeY();
        int nz0=iter->getSizeZ();

        double *dh=dhdzeta[iter-flux.begin()].getDataPtr();
        int nx1=dhdzeta[iter-flux.begin()].getSizeX();
        int ny1=dhdzeta[iter-flux.begin()].getSizeY();
        int nz1=dhdzeta[iter-flux.begin()].getSizeZ();
        
        for (int k=0; k<nz0; k++)
        {
            for (int j=0; j<ny0; j++)
            {
                for (int i=1; i<nx0; i++)
                {
                    int idx=nx0*ny0*k+nx0*j+i;
                    int idx1=nx1*ny1*k+nx1*j+i;
                    double df_loc=f[idx]-f[idx-1];
                    
                    dh[idx1-1]=df_loc;
                }
            }
        }
    }
    
}

nuc3d::EulerData3D::~EulerData3D()
{
    
}

void nuc3d::EulerData3D::solve(PDEData3d &myPDE,
                               fieldOperator3d &myOP,
                               std::vector<bufferData> &myBf,
                               physicsModel &myModel,
                               MPIComunicator3d_nonblocking &myMPI,
                               boundaryCondition &myBC)
{
    
    nuc3d::EulerData3D::solveCon2Prim(myPDE, myModel);
    nuc3d::EulerData3D::solveRiemann(myPDE, myModel);
    nuc3d::EulerData3D::setBoundaryCondition(myPDE,myModel,myBf,myBC);
    nuc3d::EulerData3D::solveInv(myOP,myBf,myMPI,myBC);
    nuc3d::EulerData3D::solveRHS(myPDE);
    
}

void nuc3d::EulerData3D::solveRiemann(PDEData3d &myPDE,
                                      physicsModel &myModel)
{
    myModel.solveRiemann(myPDE, this);
}

void nuc3d::EulerData3D::solveCon2Prim(PDEData3d &myPDE,
                                       physicsModel &myModel)
{
    myModel.solve(myPDE, this);
}

void nuc3d::EulerData3D::setBoundaryCondition(PDEData3d &myPDE,
                                              physicsModel &myModel,
                                              std::vector<bufferData> &myBf,
                                              boundaryCondition &myBC)
{
    myBC.setBC(myPDE,myModel,*this,myBf);
    
}

void nuc3d::EulerData3D::solveInv(fieldOperator3d &myOP,
                                  std::vector<bufferData> &myBf,
                                  MPIComunicator3d_nonblocking &myMPI,
                                  boundaryCondition &myBC)
{
    
    double t[8];
    
    t[0]=MPI_Wtime();
    
    solveInvicidFluxL(this->getFluxXi(), myOP, myBf,myMPI,myBC, 0);
    solveInvicidFluxL(this->getFluxEta(), myOP, myBf,myMPI,myBC, 1);
    solveInvicidFluxL(this->getFluxZeta(), myOP, myBf,myMPI,myBC, 2);
    solveInvicidFluxR(this->getFluxXi(), myOP, myBf,myMPI,myBC, 0);
    solveInvicidFluxR(this->getFluxEta(), myOP, myBf,myMPI,myBC, 1);
    solveInvicidFluxR(this->getFluxZeta(), myOP, myBf,myMPI,myBC,2);
    this->setDerivativesInv();    
}


void nuc3d::EulerData3D::solveInvicidFluxL(EulerFlux &myFlux,
                                           fieldOperator3d &myOP,
                                           std::vector<bufferData> &myBuff,
                                           MPIComunicator3d_nonblocking &myMPI,
                                           boundaryCondition &myBC,
                                           int dir)
{
    VectorField &pFlux = myFlux.FluxL;
    VectorField &pReconFlux = myFlux.reconstFluxL;
    int myBCtypeL=myBC.getBCtype(dir*2);
    int myBCtypeR=myBC.getBCtype(dir*2+1);
    
    
    for (auto iter = pFlux.begin(); iter != pFlux.end(); iter++)
    {
        bufferData &bf = myBuff[iter - pFlux.begin()];
        Field &rf = pReconFlux[iter - pFlux.begin()];
        
        myMPI.bufferSendRecv(*iter,bf,dir,static_cast<int>(iter - pFlux.begin()));
        myOP.reconstructionInner(*iter, dir, 1, rf);
        myMPI.waitSendRecv(bf,dir);
        
        Field &bfField_L=(myMPI.getFaceType(dir*2)==MPI_PROC_NULL)?bf.BufferSend[dir*2]:bf.BufferRecv[dir*2];
        Field &bfField_R=(myMPI.getFaceType(dir*2+1)==MPI_PROC_NULL)?bf.BufferSend[dir*2+1]:bf.BufferRecv[dir*2+1];
        myOP.reconstructionBoundary(*iter, bfField_L, bfField_R, dir, 1, rf,myBCtypeL,myBCtypeR);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
}

void nuc3d::EulerData3D::solveInvicidFluxR(EulerFlux &myFlux,
                                           fieldOperator3d &myOP,
                                           std::vector<bufferData> &myBuff,
                                           MPIComunicator3d_nonblocking &myMPI,
                                           boundaryCondition &myBC,
                                           int dir)

{
    VectorField &pFlux = myFlux.FluxR;
    VectorField &pReconFlux = myFlux.reconstFluxR;
    int myBCtypeL=myBC.getBCtype(dir*2);
    int myBCtypeR=myBC.getBCtype(dir*2+1);
    
    for (auto iter = pFlux.begin(); iter != pFlux.end(); iter++)
    {
        bufferData &bf = myBuff[iter - pFlux.begin()];
        Field &rf = pReconFlux[iter - pFlux.begin()];
        
        
        myMPI.bufferSendRecv(*iter,bf,dir,static_cast<int>(iter - pFlux.begin()));
        myOP.reconstructionInner(*iter, dir, -1, rf);
        
        myMPI.waitSendRecv(bf,dir);
        
        Field &bfField_L=(myMPI.getFaceType(dir*2)==MPI_PROC_NULL)?bf.BufferRecv[dir*2]:bf.BufferRecv[dir*2];
        Field &bfField_R=(myMPI.getFaceType(dir*2+1)==MPI_PROC_NULL)?bf.BufferRecv[dir*2+1]:bf.BufferRecv[dir*2+1];
        
        myOP.reconstructionBoundary(*iter, bfField_L, bfField_R, dir, -1, rf,myBCtypeL,myBCtypeR);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
}

void nuc3d::EulerData3D::solveRHS(PDEData3d &myPDE)
{
    VectorField &df_v=this->dfdxi;
    VectorField &dg_v=this->dgdeta;
    VectorField &dh_v=this->dhdzeta;
    auto beg=myPDE.getRHS().begin();
    auto end=myPDE.getRHS().end();
    
    for (auto iter=beg ; iter!=end ; iter++)
    {
        double *df=df_v[iter-beg].getDataPtr();
        double *dg=dg_v[iter-beg].getDataPtr();
        double *dh=dh_v[iter-beg].getDataPtr();
        double *rhs=iter->getDataPtr();
        
        for (int k=0; k<nz; k++)
        {
            for (int j=0; j<ny; j++)
            {
                for (int i=0; i<nx; i++)
                {
                    int idx_xi=nx*ny*k+nx*j+i;
                    
                    rhs[idx_xi]=df[idx_xi];
                }
            }
        }
        
        for (int k=0; k<nz; k++)
        {
            for (int j=0; j<ny; j++)
            {
                for (int i=0; i<nx; i++)
                {
                    int idx_xi=nx*ny*k+nx*j+i;
                    int idx_eta=ny*nz*i+ny*k+j;
                    
                    rhs[idx_xi]+=dg[idx_eta];
                }
            }
        }
        
        for (int k=0; k<nz; k++)
        {
            for (int j=0; j<ny; j++)
            {
                for (int i=0; i<nx; i++)
                {
                    int idx_xi=nx*ny*k+nx*j+i;
                    int idx_zeta=nz*nx*j+nz*i+k;
                    
                    rhs[idx_xi]+=dh[idx_zeta];
                }
            }
        }
    }
    
    getDt();
    myPDE.setDt(dt);
}

void nuc3d::EulerData3D::getDt()
{
    double dt_xi=1.0/Flux_xi.maxEigen;
    double dt_eta=1.0/Flux_eta.maxEigen;
    double dt_zeta=1.0/Flux_zeta.maxEigen;
    
    dt=std::min(std::min(dt_xi, dt_eta),dt_zeta);}

#endif