#ifndef euler3d_cpp
#define euler3d_cpp
#include <cmath>
#include "euler3d.h"
#include "physicsModel.h"
#include "PDEData3d.hpp"

/**************************************************************************************
 Member functions of class: EulerFlux
 **************************************************************************************/
nuc3d::EulerFlux::EulerFlux(int nx0,int ny0,int nz0,int neqs,
                            int xdir,int ydir,int zdir):
FluxL(neqs,Field(nx0,ny0,nz0)),
FluxR(neqs,Field(nx0,ny0,nz0)),
reconstFluxL(neqs,Field(nx0+xdir,ny0+ydir,nz0+zdir)),
reconstFluxR(neqs,Field(nx0+xdir,ny0+ydir,nz0+zdir)),
reconstFlux(neqs,Field(nx0,ny0,nz0)),
maxEigen(1.0)
{
    
}

void nuc3d::EulerFlux::combineFluxLR()
{
    for(auto iter=reconstFlux.begin();iter!=reconstFlux.end();iter++)
        *iter=reconstFluxL[iter-reconstFlux.begin()]
        +reconstFluxR[iter-reconstFlux.begin()];
}

nuc3d::EulerFlux::~EulerFlux()
{}

/**************************************************************************************
 Member functions of class: EulerData3D
 **************************************************************************************/
nuc3d::EulerData3D::EulerData3D( int nx0, int ny0, int nz0, int neqs):
jacobian(Field(nx0,ny0,nz0)),
xi_xyz(3,Field(nx0,ny0,nz0)),
eta_xyz(3,Field(nx0,ny0,nz0)),
zeta_xyz(3,Field(nx0,ny0,nz0)),
W0_Euler(3,Field(nx0,ny0,nz0)),
W_Euler(neqs,Field(nx0,ny0,nz0)),
Flux_xi(nx0,ny0,nz0,neqs,1,0,0),
Flux_eta(nx0,ny0,nz0,neqs,0,1,0),
Flux_zeta(nx0,ny0,nz0,neqs,0,0,1),
dfdxi(neqs,Field(nx0,ny0,nz0)),
dgdeta(neqs,Field(nx0,ny0,nz0)),
dhdzeta(neqs,Field(nx0,ny0,nz0)),
dt(0.0)
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
    
    setDerivativesXiInv();
    setDerivativesEtaInv();
    setDerivativesZetaInv();
}

void nuc3d::EulerData3D::setDerivativesXiInv()
{
    VectorField &flux=Flux_xi.reconstFlux;
    for (auto iter=flux.begin() ; iter!=flux.end(); iter++)
    {
        int nx0=iter->getSizeX();
        int ny0=iter->getSizeY();
        int nz0=iter->getSizeZ();
        
        for (int k=0; k<nz0; k++)
        {
            for (int j=0; j<ny0; j++)
            {
                for (int i=1; i<nx0; i++)
                {
                    double df=iter->getValue(i, j, k)-iter->getValue(i-1, j, k);
                    
                    dfdxi[iter-flux.begin()].setValue(i-1, j, k, df);
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
        int nx0=iter->getSizeX();
        int ny0=iter->getSizeY();
        int nz0=iter->getSizeZ();
        
        for (int k=0; k<nz0; k++)
        {
            for (int j=1; j<ny0; j++)
            {
                for (int i=0; i<nx0; i++)
                {
                    double df=iter->getValue(i, j, k)-iter->getValue(i, j-1, k);
                    
                    dgdeta[iter-flux.begin()].setValue(i, j-1, k, df);
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
        int nx0=iter->getSizeX();
        int ny0=iter->getSizeY();
        int nz0=iter->getSizeZ();
        
        for (int k=1; k<nz0; k++)
        {
            for (int j=0; j<ny0; j++)
            {
                for (int i=0; i<nx0; i++)
                {
                    double df=iter->getValue(i, j, k)-iter->getValue(i, j, k-1);
                    
                    dfdxi[iter-flux.begin()].setValue(i, j, k-1, df);
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
                               MPIComunicator3d_nonblocking &myMPI)
{
    nuc3d::EulerData3D::solveRiemann(myPDE, myModel);
    nuc3d::EulerData3D::solveInv(myOP,myBf,myMPI);
    nuc3d::EulerData3D::solveRHS(myPDE);
}

void nuc3d::EulerData3D::solveRiemann(PDEData3d &myPDE,
                                      physicsModel &myModel)
{
    myModel.solve(myPDE, this);
}

void nuc3d::EulerData3D::solveInv(fieldOperator3d &myOP,
                                  std::vector<bufferData> &myBf,
                                  MPIComunicator3d_nonblocking &myMPI)
{
    solveInvicidFluxL(this->getFluxXi(), myOP, myBf,myMPI, 0);
    solveInvicidFluxL(this->getFluxEta(), myOP, myBf,myMPI, 1);
    solveInvicidFluxL(this->getFluxZeta(), myOP, myBf,myMPI, 2);
    
    solveInvicidFluxR(this->getFluxXi(), myOP, myBf,myMPI, 0);
    solveInvicidFluxR(this->getFluxEta(), myOP, myBf,myMPI, 1);
    solveInvicidFluxR(this->getFluxZeta(), myOP, myBf,myMPI, 2);
    
    this->setDerivativesInv();
}


void nuc3d::EulerData3D::solveInvicidFluxL(EulerFlux &myFlux,
                                           fieldOperator3d &myOP,
                                           std::vector<bufferData> &myBuff,
                                           MPIComunicator3d_nonblocking &myMPI,
                                           int dir)
{
    VectorField &pFlux = myFlux.FluxL;
    VectorField &pReconFlux = myFlux.reconstFluxL;
    
    for (auto iter = pFlux.begin(); iter != pFlux.end(); iter++)
    {
        Field &rf = pReconFlux[iter - pFlux.begin()];
        bufferData &bf = myBuff[iter - pFlux.begin()];
        
        bf.setBufferSend(*iter);
        
        myMPI.bufferSendRecv(bf, 0);
        
        myOP.reconstructionInner(*iter, dir, 1, rf);
        
        myMPI.waitAllSendRecv(bf);
        
        myOP.reconstructionBoundary(*iter, bf.BufferRecv[dir], bf.BufferRecv[dir+1], dir, 1, rf);
        
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
}

void nuc3d::EulerData3D::solveInvicidFluxR(EulerFlux &myFlux,
                                           fieldOperator3d &myOP,
                                           std::vector<bufferData> &myBuff,
                                           MPIComunicator3d_nonblocking &myMPI,
                                           int dir)

{
    VectorField &pFlux = myFlux.FluxR;
    VectorField &pReconFlux = myFlux.reconstFluxR;
    
    for (auto iter = pFlux.begin(); iter != pFlux.end(); iter++)
    {
        Field &rf = pReconFlux[iter - pFlux.begin()];
        bufferData &bf = myBuff[iter - pFlux.begin()];
        
        bf.setBufferSend(*iter);
        
        myMPI.bufferSendRecv(bf, 0);
        
        myOP.reconstructionInner(*iter, dir, -1, rf);
        
        myMPI.waitAllSendRecv(bf);
        
        myOP.reconstructionBoundary(*iter, bf.BufferRecv[dir], bf.BufferRecv[dir + 1], dir, -1, rf);
        
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
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
        int nx0=iter->getSizeX();
        int ny0=iter->getSizeY();
        int nz0=iter->getSizeZ();
        
        for (int k=0; k<nz0; k++)
        {
            for (int j=0; j<ny0; j++)
            {
                for (int i=0; i<nx0; i++)
                {
                    double df=df_v[iter-beg].getValue(i, j, k);
                    double dg=dg_v[iter-beg].getValue(i, j, k);
                    double dh=dh_v[iter-beg].getValue(i, j, k);
                    
                    double rhs=df+dg+dh;
                    
                    iter->setValue(i, j, k, rhs);
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
    
    dt=std::min(std::min(dt_xi, dt_eta),dt_zeta);
}

#endif