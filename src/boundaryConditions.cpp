//
//  boundaryConditions.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/11/5.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#include "boundaryConditions.hpp"
#include "MPICommunicator.h"
#include "euler3d.h"
#include "NaiverStokes3d.h"
#include "PDEData3d.hpp"
#include "physicsModel.h"

nuc3d::boundaryCondition::boundaryCondition()
{
}

nuc3d::boundaryCondition::~boundaryCondition()
{
    
}

void nuc3d::boundaryCondition::initialBC(VectorBuffer &myBuffer,
                                         MPIComunicator3d_nonblocking &myMPI)
{
    std::stringstream ss;
    int ProcId = myMPI.getMyId();
    std::string forename_bc = ("bc_topo/BC_Topo_");
    std::string midname;
    std::string tailname = (".in");
    
    ss << ProcId;
    ss >> midname;
    
    std::string filename_bc = forename_bc + midname + tailname;
    
    std::ifstream myFile;
    myFile.open(filename_bc);
    
    if (myFile)
    {
        int nline=0;
        while(readBCTopo(myFile,nline)&&(6!=nline))
        {
            nline++;
        };
        
        if(6!=BCTopo.size())
        {
            std::cout<<"Error: BC file:"<<filename_bc<<" is incomplete"<<std::endl;
            std::cout<<"       BC number is less than 6:"<<BCTopo.size()<<std::endl;
            exit(-1);
        }
    }
    else
    {
        std::cout<<"file "<<filename_bc<<" does not exist!"<<std::endl;
        exit(-1);
    }
    
    myFile.close();
    myMPI.setTopo(BCTopo);
    
    
}


std::ifstream& nuc3d::boundaryCondition::readBCTopo(std::ifstream& ios,int iface)
{
    std::string line;
    std::vector<double> BCvalue0(0);
    if(std::getline(ios,line))
    {
        int type;
        int id;
        double value;
        std::istringstream values(line);
        values>>type>>id;
        BCTopo.push_back(faceBC(type,id));
        
        while(values>>value)
        {
            BCvalue0.push_back(value);
        }
        
        
        BCvalue.push_back(BCvalue0);
        
    }
    
    return ios;
}

void nuc3d::boundaryCondition::setBC(PDEData3d &myPDE,
                                     physicsModel &myPhyMod,
                                     EulerData3D &myFluxes,
                                     VectorBuffer &myBf)
{
    for(auto iter=BCTopo.begin();iter!=BCTopo.end();iter++)
    {
        int iface=static_cast<int>(iter-BCTopo.begin());
        if (((-1)==iter->Type))
        {
            (this->*mySetter[BCTopo[iface].id])(myPDE,myPhyMod,myFluxes,myBf,iface);
        }
    }
}
//Inlet condition
void nuc3d::boundaryCondition::setBC_Inlet(PDEData3d &myPDE,
                                           physicsModel &myPhyMod,
                                           EulerData3D &myFluxes,
                                           VectorBuffer &myBf,
                                           int iface)
{
    switch (iface) {
        case 0:
            BCsetter_inlet_xi(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 1:
            BCsetter_inlet_xi(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 2:
            BCsetter_inlet_eta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 3:
            BCsetter_inlet_eta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 4:
            BCsetter_inlet_zeta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 5:
            BCsetter_inlet_zeta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
    }
    
}

void nuc3d::boundaryCondition::BCsetter_inlet_xi(PDEData3d &myPDE,
                                                 physicsModel &myPhyMod,
                                                 EulerData3D &myFluxes,
                                                 VectorBuffer &myBf,
                                                 int lr)
{
    std::vector<double> fluxl(myPhyMod.getEqNum());
    std::vector<double> fluxr(myPhyMod.getEqNum());
    
    Field &jac=myFluxes.getJac();
    VectorField &xi_xyz=myFluxes.getXi_xyz();
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(int k=0;k<nz;k++)
    {
        for(int j=0;j<ny;j++)
        {
            double jacob=jac.getValue(lr*(nx-1),j,k);
            double xi_x=xi_xyz[0].getValue(lr*(nx-1),j,k);
            double xi_y=xi_xyz[1].getValue(lr*(nx-1),j,k);
            double xi_z=xi_xyz[2].getValue(lr*(nx-1),j,k);
            
            myPhyMod.solveRiemannPoint(BCvalue[lr], jacob, xi_x, xi_y, xi_z, fluxl, fluxr);
            
            for(auto iter=myBf.begin();iter!=myBf.end();iter++)
            {
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    
                    iter->BufferSend[lr].setValue(ibf, j, k, fluxl[iter-myBf.begin()]);
                    iter->BufferRecv[lr].setValue(ibf, j, k, fluxr[iter-myBf.begin()]);
                }
            }
        }
    }
    
}

void nuc3d::boundaryCondition::BCsetter_inlet_eta(PDEData3d &myPDE,
                                                  physicsModel &myPhyMod,
                                                  EulerData3D &myFluxes,
                                                  VectorBuffer &myBf,
                                                  int lr)
{
    std::vector<double> fluxl(myPhyMod.getEqNum());
    std::vector<double> fluxr(myPhyMod.getEqNum());
    
    Field &jac=myFluxes.getJac();
    VectorField &eta_xyz=myFluxes.getEta_xyz();
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(int k=0;k<nz;k++)
    {
        for(int i=0;i<nx;i++)
        {
            double jacob=jac.getValue(i,lr*(ny-1),k);
            double eta_x=eta_xyz[0].getValue(i,lr*(ny-1),k);
            double eta_y=eta_xyz[1].getValue(i,lr*(ny-1),k);
            double eta_z=eta_xyz[2].getValue(i,lr*(ny-1),k);
            
            myPhyMod.solveRiemannPoint(BCvalue[2+lr], jacob, eta_x, eta_y, eta_z, fluxl, fluxr);
            
            for(auto iter=myBf.begin();iter!=myBf.end();iter++)
            {
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    
                    iter->BufferSend[2+lr].setValue(i,ibf, k, fluxl[iter-myBf.begin()]);
                    iter->BufferRecv[2+lr].setValue(i,ibf, k, fluxr[iter-myBf.begin()]);
                }
            }
        }
    }
}

void nuc3d::boundaryCondition::BCsetter_inlet_zeta(PDEData3d &myPDE,
                                                   physicsModel &myPhyMod,
                                                   EulerData3D &myFluxes,
                                                   VectorBuffer &myBf,
                                                   int lr)
{
    std::vector<double> fluxl(myPhyMod.getEqNum());
    std::vector<double> fluxr(myPhyMod.getEqNum());
    
    Field &jac=myFluxes.getJac();
    VectorField &zeta_xyz=myFluxes.getZeta_xyz();
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(int j=0;j<ny;j++)
    {
        for(int i=0;i<nx;i++)
        {
            double jacob=jac.getValue(i,j,lr*(nz-1));
            double zeta_x=zeta_xyz[0].getValue(i,j,lr*(nz-1));
            double zeta_y=zeta_xyz[1].getValue(i,j,lr*(nz-1));
            double zeta_z=zeta_xyz[2].getValue(i,j,lr*(nz-1));
            
            myPhyMod.solveRiemannPoint(BCvalue[4+lr], jacob, zeta_x, zeta_y, zeta_z, fluxl, fluxr);
            
            for(auto iter=myBf.begin();iter!=myBf.end();iter++)
            {
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    
                    iter->BufferSend[4+lr].setValue(i,j,ibf, fluxl[iter-myBf.begin()]);
                    iter->BufferRecv[4+lr].setValue(i,j,ibf, fluxr[iter-myBf.begin()]);
                }
            }
        }
    }
}

// Outlet condition
void nuc3d::boundaryCondition::setBC_Outlet(PDEData3d &myPDE,
                                            physicsModel &myPhyMod,
                                            EulerData3D &myFluxes,
                                            VectorBuffer &myBf,
                                            int iface)
{
    switch (iface) {
        case 0:
            BCsetter_outlet_xi(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 1:
            BCsetter_outlet_xi(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 2:
            BCsetter_outlet_eta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 3:
            BCsetter_outlet_eta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 4:
            BCsetter_outlet_zeta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 5:
            BCsetter_outlet_zeta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
    }
    
    
}

void nuc3d::boundaryCondition::BCsetter_outlet_xi(PDEData3d &myPDE,
                                                  physicsModel &myPhyMod,
                                                  EulerData3D &myFluxes,
                                                  VectorBuffer &myBf,
                                                  int lr)
{
    std::vector<double> fluxl(myPhyMod.getEqNum());
    std::vector<double> fluxr(myPhyMod.getEqNum());
    std::vector<double> q(myPhyMod.getEqNum());
    
    
    Field &jac=myFluxes.getJac();
    VectorField &xi_xyz=myFluxes.getXi_xyz();
    VectorField &prim=myFluxes.getPrimatives();
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(int k=0;k<nz;k++)
    {
        for(int j=0;j<ny;j++)
        {
            double jacob=jac.getValue(lr*(nx-1),j,k);
            double xi_x=xi_xyz[0].getValue(lr*(nx-1),j,k);
            double xi_y=xi_xyz[1].getValue(lr*(nx-1),j,k);
            double xi_z=xi_xyz[2].getValue(lr*(nx-1),j,k);
            
            for(auto iter=q.begin();iter!=q.end();iter++)
                *iter=prim[iter-q.begin()].getValue(lr*(nx-1), j, k);
            
            myPhyMod.solveRiemannPoint(q, jacob, xi_x, xi_y, xi_z, fluxl, fluxr);
            
            for(auto iter=myBf.begin();iter!=myBf.end();iter++)
            {
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    
                    iter->BufferSend[lr].setValue(ibf, j, k, fluxl[iter-myBf.begin()]);
                    iter->BufferRecv[lr].setValue(ibf, j, k, fluxr[iter-myBf.begin()]);
                }
            }
        }
    }
    
}

void nuc3d::boundaryCondition::BCsetter_outlet_eta(PDEData3d &myPDE,
                                                   physicsModel &myPhyMod,
                                                   EulerData3D &myFluxes,
                                                   VectorBuffer &myBf,
                                                   int lr)
{
    std::vector<double> fluxl(myPhyMod.getEqNum());
    std::vector<double> fluxr(myPhyMod.getEqNum());
    std::vector<double> q(myPhyMod.getEqNum());
    
    Field &jac=myFluxes.getJac();
    VectorField &eta_xyz=myFluxes.getEta_xyz();
    VectorField &prim=myFluxes.getPrimatives();
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(int k=0;k<nz;k++)
    {
        for(int i=0;i<nx;i++)
        {
            double jacob=jac.getValue(i,lr*(ny-1),k);
            double eta_x=eta_xyz[0].getValue(i,lr*(ny-1),k);
            double eta_y=eta_xyz[1].getValue(i,lr*(ny-1),k);
            double eta_z=eta_xyz[2].getValue(i,lr*(ny-1),k);
            
            for(auto iter=q.begin();iter!=q.end();iter++)
                *iter=prim[iter-q.begin()].getValue(i,lr*(ny-1), k);
            
            myPhyMod.solveRiemannPoint(q, jacob, eta_x, eta_y, eta_z, fluxl, fluxr);
            
            for(auto iter=myBf.begin();iter!=myBf.end();iter++)
            {
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    
                    iter->BufferSend[2+lr].setValue(i,ibf, k, fluxl[iter-myBf.begin()]);
                    iter->BufferRecv[2+lr].setValue(i,ibf, k, fluxr[iter-myBf.begin()]);
                }
            }
        }
    }
    
}
void nuc3d::boundaryCondition::BCsetter_outlet_zeta(PDEData3d &myPDE,
                                                    physicsModel &myPhyMod,
                                                    EulerData3D &myFluxes,
                                                    VectorBuffer &myBf,
                                                    int lr)
{
    std::vector<double> fluxl(myPhyMod.getEqNum());
    std::vector<double> fluxr(myPhyMod.getEqNum());
    std::vector<double> q(myPhyMod.getEqNum());
    
    Field &jac=myFluxes.getJac();
    VectorField &zeta_xyz=myFluxes.getZeta_xyz();
    VectorField &prim=myFluxes.getPrimatives();
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(int k=0;k<bfsize;k++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int i=0;i<nx;i++)
            {
                int kbuff=(bfsize-1)*(1-lr)-(1-2*lr)*k;
                
                double jacob=jac.getValue(i,j,lr*(nz-1));
                
                double zeta_x=zeta_xyz[0].getValue(i,j,lr*(nz-1));
                double zeta_y=zeta_xyz[1].getValue(i,j,lr*(nz-1));
                double zeta_z=zeta_xyz[2].getValue(i,j,lr*(nz-1));
                
                for(auto iter=q.begin();iter!=q.end();iter++)
                    *iter=prim[iter-q.begin()].getValue(i,j,lr*(nz-1));
                
                myPhyMod.solveRiemannPoint(q, jacob, zeta_x, zeta_y, zeta_z, fluxl, fluxr);
                
                for(auto iter=myBf.begin();iter!=myBf.end();iter++)
                {
                    iter->BufferSend[4+lr].setValue(i,j,kbuff, fluxl[iter-myBf.begin()]);
                    iter->BufferRecv[4+lr].setValue(i,j,kbuff, fluxr[iter-myBf.begin()]);
                }
            }
            
        }
    }
    
}

void nuc3d::boundaryCondition::setBC_wall(PDEData3d &myPDE,
                                          physicsModel &myPhyMod,
                                          EulerData3D &myFluxes,
                                          VectorBuffer &myBf,
                                          int iface)
{
    switch (iface) {
        case 0:
            BCsetter_wall_xi(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 1:
            BCsetter_wall_xi(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 2:
            BCsetter_wall_eta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 3:
            BCsetter_wall_eta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 4:
            BCsetter_wall_zeta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 5:
            BCsetter_wall_zeta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
    }
}

void nuc3d::boundaryCondition::BCsetter_wall_xi(PDEData3d &myPDE,
                                                physicsModel &myPhyMod,
                                                EulerData3D &myFluxes,
                                                VectorBuffer &myBf,
                                                int lr)
{
    std::vector<double> fluxl(myPhyMod.getEqNum());
    std::vector<double> fluxr(myPhyMod.getEqNum());
    std::vector<double> q(myPhyMod.getEqNum());
    
    Field &jac=myFluxes.getJac();
    VectorField &xi_xyz=myFluxes.getXi_xyz();
    VectorField &prim=myFluxes.getPrimatives();
    VectorField &accu=myFluxes.getAcoustics();
    double T;
    
    double T_wall=myPhyMod.getWallTemp();
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(int k=0;k<nz;k++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int i=0;i<bfsize;i++)
            {
                int ibuff=(bfsize-1)*(1-lr)-(1-2*lr)*i;
                
                double jacob=jac.getValue(lr*(nx-1),j,k);
                double xi_x=xi_xyz[0].getValue(lr*(nx-1),j,k);
                double xi_y=xi_xyz[1].getValue(lr*(nx-1),j,k);
                double xi_z=xi_xyz[2].getValue(lr*(nx-1),j,k);
                
                for(auto iter=q.begin();iter!=q.end();iter++)
                    *iter=prim[iter-q.begin()].getValue(lr*(nx-1), j, k);
                
                T=accu[0].getValue(lr*(nx-1), j, k);
                
                q[0]=q[0]*T/T_wall;
                q[1]=0.0;
                q[2]=0.0;
                q[3]=0.0;
                
                myPhyMod.solveRiemannPoint(q, jacob, xi_x, xi_y, xi_z, fluxl, fluxr);
                
                for(auto iter=myBf.begin();iter!=myBf.end();iter++)
                {
                    iter->BufferSend[lr].setValue(ibuff, j, k, fluxl[iter-myBf.begin()]);
                    iter->BufferRecv[lr].setValue(ibuff, j, k, fluxr[iter-myBf.begin()]);
                }
            }
            
        }
    }
}

void nuc3d::boundaryCondition::BCsetter_wall_eta(PDEData3d &myPDE,
                                                 physicsModel &myPhyMod,
                                                 EulerData3D &myFluxes,
                                                 VectorBuffer &myBf,
                                                 int lr)
{
    std::vector<double> fluxl(myPhyMod.getEqNum());
    std::vector<double> fluxr(myPhyMod.getEqNum());
    std::vector<double> q(myPhyMod.getEqNum());
    Field &jac=myFluxes.getJac();
    VectorField &eta_xyz=myFluxes.getEta_xyz();
    
    VectorField &prim=myFluxes.getPrimatives();
    VectorField &accu=myFluxes.getAcoustics();
    double T;
    double T_wall=myPhyMod.getWallTemp();
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(int k=0;k<nz;k++)
    {
        for(int j=0;j<bfsize;j++)
        {
            for(int i=0;i<nx;i++)
            {
                int jbuff=(bfsize-1)*(1-lr)-(1-2*lr)*j;
                
                double jacob=jac.getValue(i,lr*(ny-1),k);
                
                double eta_x=eta_xyz[0].getValue(i,lr*(ny-1),k);
                double eta_y=eta_xyz[1].getValue(i,lr*(ny-1),k);
                double eta_z=eta_xyz[2].getValue(i,lr*(ny-1),k);
                
                
                for(auto iter=q.begin();iter!=q.end();iter++)
                    *iter=prim[iter-q.begin()].getValue(i,lr*(ny-1), k);
                
                T=accu[0].getValue(i,lr*(ny-1), k);
                q[0]=q[0]*T/T_wall;
                q[1]=0.0;
                q[2]=0.0;
                q[3]=0.0;
                
                myPhyMod.solveRiemannPoint(q, jacob, eta_x, eta_y, eta_z, fluxl, fluxr);
                
                for(auto iter=myBf.begin();iter!=myBf.end();iter++)
                {
                    iter->BufferSend[2+lr].setValue(i,jbuff, k, fluxl[iter-myBf.begin()]);
                    iter->BufferRecv[2+lr].setValue(i,jbuff, k, fluxr[iter-myBf.begin()]);
                }
            }
            
        }
    }
}

void nuc3d::boundaryCondition::BCsetter_wall_zeta(PDEData3d &myPDE,
                                                  physicsModel &myPhyMod,
                                                  EulerData3D &myFluxes,
                                                  VectorBuffer &myBf,
                                                  int lr)
{
    std::vector<double> fluxl(myPhyMod.getEqNum());
    std::vector<double> fluxr(myPhyMod.getEqNum());
    std::vector<double> q(myPhyMod.getEqNum());
    
    Field &jac=myFluxes.getJac();
    VectorField &zeta_xyz=myFluxes.getZeta_xyz();
    
    VectorField &prim=myFluxes.getPrimatives();
    VectorField &accu=myFluxes.getAcoustics();
    double T;
    double T_wall=myPhyMod.getWallTemp();
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(int k=0;k<bfsize;k++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int i=0;i<nx;i++)
            {
                int kbuff=(bfsize-1)*(1-lr)-(1-2*lr)*k;
                
                double jacob=jac.getValue(i,j,lr*(nz-1));
                
                double zeta_x=zeta_xyz[0].getValue(i,j,lr*(nz-1));
                double zeta_y=zeta_xyz[1].getValue(i,j,lr*(nz-1));
                double zeta_z=zeta_xyz[2].getValue(i,j,lr*(nz-1));
                
                
                for(auto iter=q.begin();iter!=q.end();iter++)
                    *iter=prim[iter-q.begin()].getValue(i,j,lr*(nz-1));
                
                T=accu[0].getValue(i,j,lr*(nz-1));
                
                q[0]=q[0]*T/T_wall;
                
                q[1]=0.0;
                q[2]=0.0;
                q[3]=0.0;
                
                myPhyMod.solveRiemannPoint(q, jacob, zeta_x, zeta_y, zeta_z, fluxl, fluxr);
                
                for(auto iter=myBf.begin();iter!=myBf.end();iter++)
                {
                    iter->BufferSend[4+lr].setValue(i,j,kbuff, fluxl[iter-myBf.begin()]);
                    iter->BufferRecv[4+lr].setValue(i,j,kbuff, fluxr[iter-myBf.begin()]);
                }
            }
            
        }
    }
}


void nuc3d::boundaryCondition::setBC_symm(PDEData3d &myPDE,
                                          physicsModel &myPhyMod,
                                          EulerData3D &myFluxes,
                                          VectorBuffer &myBf,
                                          int iface)
{
    switch (iface) {
        case 0:
            BCsetter_symm_xi(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 1:
            BCsetter_symm_xi(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 2:
            BCsetter_symm_eta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 3:
            BCsetter_symm_eta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 4:
            BCsetter_symm_zeta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 5:
            BCsetter_symm_zeta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
    }
    
}

void nuc3d::boundaryCondition::BCsetter_symm_xi(PDEData3d &myPDE,
                                                physicsModel &myPhyMod,
                                                EulerData3D &myFluxes,
                                                VectorBuffer &myBf,
                                                int lr)
{
    std::vector<double> fluxl(myPhyMod.getEqNum());
    std::vector<double> fluxr(myPhyMod.getEqNum());
    std::vector<double> q(myPhyMod.getEqNum());
    
    Field &jac=myFluxes.getJac();
    VectorField &xi_xyz=myFluxes.getXi_xyz();
    VectorField &eta_xyz=myFluxes.getEta_xyz();
    VectorField &zeta_xyz=myFluxes.getZeta_xyz();
    
    VectorField &dx=myFluxes.getDx();
    VectorField &dy=myFluxes.getDy();
    VectorField &dz=myFluxes.getDz();
    
    VectorField &prim=myFluxes.getPrimatives();
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(int k=0;k<nz;k++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int i=0;i<bfsize;i++)
            {
                int iblock=lr*(nx-1)+(1-2*lr)*i;
                int ibuff=(bfsize-1)*(1-lr)-(1-2*lr)*i;
                
                double jacob=jac.getValue(iblock,j,k);
                double xi_x=xi_xyz[0].getValue(iblock,j,k);
                double xi_y=xi_xyz[1].getValue(iblock,j,k);
                double xi_z=xi_xyz[2].getValue(iblock,j,k);
                
                double eta_x=eta_xyz[0].getValue(iblock,j,k);
                double eta_y=eta_xyz[1].getValue(iblock,j,k);
                double eta_z=eta_xyz[2].getValue(iblock,j,k);
                
                double zeta_x=zeta_xyz[0].getValue(iblock,j,k);
                double zeta_y=zeta_xyz[1].getValue(iblock,j,k);
                double zeta_z=zeta_xyz[2].getValue(iblock,j,k);
                
                double x_xi=dx[0].getValue(iblock,j,k);
                double y_xi=dy[0].getValue(iblock,j,k);
                double z_xi=dz[0].getValue(iblock,j,k);
                
                double x_eta=dx[1].getValue(iblock,j,k);
                double y_eta=dy[1].getValue(iblock,j,k);
                double z_eta=dz[1].getValue(iblock,j,k);
                
                double x_zeta=dx[2].getValue(iblock,j,k);
                double y_zeta=dy[2].getValue(iblock,j,k);
                double z_zeta=dz[2].getValue(iblock,j,k);
                
                for(auto iter=q.begin();iter!=q.end();iter++)
                    *iter=prim[iter-q.begin()].getValue(iblock, j, k);
                
                double U=xi_x*q[1]+xi_y*q[2]+xi_z*q[3];
                double V=eta_x*q[1]+eta_y*q[2]+eta_z*q[3];
                double W=zeta_x*q[1]+zeta_y*q[2]+zeta_z*q[3];
                
                q[1]=-x_xi*U+x_eta*V+x_zeta*W;
                q[2]=-y_xi*U+y_eta*V+y_zeta*W;
                q[3]=-z_xi*U+z_eta*V+z_zeta*W;
                
                myPhyMod.solveRiemannPoint(q, jacob, xi_x, xi_y, xi_z, fluxl, fluxr);
                
                for(auto iter=myBf.begin();iter!=myBf.end();iter++)
                {
                    iter->BufferSend[lr].setValue(ibuff, j, k, fluxl[iter-myBf.begin()]);
                    iter->BufferRecv[lr].setValue(ibuff, j, k, fluxr[iter-myBf.begin()]);
                }
            }
            
        }
    }
}

void nuc3d::boundaryCondition::BCsetter_symm_eta(PDEData3d &myPDE,
                                                 physicsModel &myPhyMod,
                                                 EulerData3D &myFluxes,
                                                 VectorBuffer &myBf,
                                                 int lr)
{
    std::vector<double> fluxl(myPhyMod.getEqNum());
    std::vector<double> fluxr(myPhyMod.getEqNum());
    std::vector<double> q(myPhyMod.getEqNum());
    
    Field &jac=myFluxes.getJac();
    VectorField &xi_xyz=myFluxes.getXi_xyz();
    VectorField &eta_xyz=myFluxes.getEta_xyz();
    VectorField &zeta_xyz=myFluxes.getZeta_xyz();
    
    VectorField &dx=myFluxes.getDx();
    VectorField &dy=myFluxes.getDy();
    VectorField &dz=myFluxes.getDz();
    
    VectorField &prim=myFluxes.getPrimatives();
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(int k=0;k<nz;k++)
    {
        for(int j=0;j<bfsize;j++)
        {
            for(int i=0;i<nx;i++)
            {
                int jblock=lr*(ny-1)+(1-2*lr)*j;
                int jbuff=(bfsize-1)*(1-lr)-(1-2*lr)*j;
                
                double jacob=jac.getValue(i,jblock,k);
                
                double xi_x=xi_xyz[0].getValue(i,jblock,k);
                double xi_y=xi_xyz[1].getValue(i,jblock,k);
                double xi_z=xi_xyz[2].getValue(i,jblock,k);
                
                double eta_x=eta_xyz[0].getValue(i,jblock,k);
                double eta_y=eta_xyz[1].getValue(i,jblock,k);
                double eta_z=eta_xyz[2].getValue(i,jblock,k);
                
                double zeta_x=zeta_xyz[0].getValue(i,jblock,k);
                double zeta_y=zeta_xyz[1].getValue(i,jblock,k);
                double zeta_z=zeta_xyz[2].getValue(i,jblock,k);
                
                double x_xi=dx[0].getValue(i, jblock, k);
                double y_xi=dy[0].getValue(i, jblock, k);
                double z_xi=dz[0].getValue(i, jblock, k);
                
                double x_eta=dx[1].getValue(i, jblock, k);
                double y_eta=dy[1].getValue(i, jblock, k);
                double z_eta=dz[1].getValue(i, jblock, k);
                
                double x_zeta=dx[2].getValue(i, jblock, k);
                double y_zeta=dy[2].getValue(i, jblock, k);
                double z_zeta=dz[2].getValue(i, jblock, k);
                
                for(auto iter=q.begin();iter!=q.end();iter++)
                    *iter=prim[iter-q.begin()].getValue(i,jblock, k);
                
                double U=xi_x*q[1]+xi_y*q[2]+xi_z*q[3];
                double V=eta_x*q[1]+eta_y*q[2]+eta_z*q[3];
                double W=zeta_x*q[1]+zeta_y*q[2]+zeta_z*q[3];
                q[1]=x_xi*U-x_eta*V+x_zeta*W;
                q[2]=y_xi*U-y_eta*V+y_zeta*W;
                q[3]=z_xi*U-z_eta*V+z_zeta*W;
                
                
                myPhyMod.solveRiemannPoint(q, jacob, eta_x, eta_y, eta_z, fluxl, fluxr);
                
                for(auto iter=myBf.begin();iter!=myBf.end();iter++)
                {
                    iter->BufferSend[2+lr].setValue(i,jbuff, k, fluxl[iter-myBf.begin()]);
                    iter->BufferRecv[2+lr].setValue(i,jbuff, k, fluxr[iter-myBf.begin()]);
                }
            }
            
        }
    }
}

void nuc3d::boundaryCondition::BCsetter_symm_zeta(PDEData3d &myPDE,
                                                  physicsModel &myPhyMod,
                                                  EulerData3D &myFluxes,
                                                  VectorBuffer &myBf,
                                                  int lr)
{
    std::vector<double> fluxl(myPhyMod.getEqNum());
    std::vector<double> fluxr(myPhyMod.getEqNum());
    std::vector<double> q(myPhyMod.getEqNum());
    
    Field &jac=myFluxes.getJac();
    VectorField &xi_xyz=myFluxes.getXi_xyz();
    VectorField &eta_xyz=myFluxes.getEta_xyz();
    VectorField &zeta_xyz=myFluxes.getZeta_xyz();
    
    VectorField &dx=myFluxes.getDx();
    VectorField &dy=myFluxes.getDy();
    VectorField &dz=myFluxes.getDz();
    
    VectorField &prim=myFluxes.getPrimatives();
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(int k=0;k<bfsize;k++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int i=0;i<nx;i++)
            {
                int kblock=lr*(nz-1)+(1-2*lr)*k;
                int kbuff=(bfsize-1)*(1-lr)-(1-2*lr)*k;
                
                double jacob=jac.getValue(i,j,kblock);
                
                double xi_x=xi_xyz[0].getValue(i,j,kblock);
                double xi_y=xi_xyz[1].getValue(i,j,kblock);
                double xi_z=xi_xyz[2].getValue(i,j,kblock);
                
                double eta_x=eta_xyz[0].getValue(i,j,kblock);
                double eta_y=eta_xyz[1].getValue(i,j,kblock);
                double eta_z=eta_xyz[2].getValue(i,j,kblock);
                
                double zeta_x=zeta_xyz[0].getValue(i,j,kblock);
                double zeta_y=zeta_xyz[1].getValue(i,j,kblock);
                double zeta_z=zeta_xyz[2].getValue(i,j,kblock);
                
                double x_xi=dx[0].getValue(i,j,kblock);
                double y_xi=dy[0].getValue(i,j,kblock);
                double z_xi=dz[0].getValue(i,j,kblock);
                
                double x_eta=dx[1].getValue(i,j,kblock);
                double y_eta=dy[1].getValue(i,j,kblock);
                double z_eta=dz[1].getValue(i,j,kblock);
                
                double x_zeta=dx[2].getValue(i,j,kblock);
                double y_zeta=dy[2].getValue(i,j,kblock);
                double z_zeta=dz[2].getValue(i,j,kblock);
                
                for(auto iter=q.begin();iter!=q.end();iter++)
                    *iter=prim[iter-q.begin()].getValue(i,j,kblock);
                
                double U=xi_x*q[1]+xi_y*q[2]+xi_z*q[3];
                double V=eta_x*q[1]+eta_y*q[2]+eta_z*q[3];
                double W=zeta_x*q[1]+zeta_y*q[2]+zeta_z*q[3];
                
                q[1]=x_xi*U+x_eta*V-x_zeta*W;
                q[2]=y_xi*U+y_eta*V-y_zeta*W;
                q[3]=z_xi*U+z_eta*V-z_zeta*W;
                
                
                myPhyMod.solveRiemannPoint(q, jacob, zeta_x, zeta_y, zeta_z, fluxl, fluxr);
                
                for(auto iter=myBf.begin();iter!=myBf.end();iter++)
                {
                    iter->BufferSend[4+lr].setValue(i,j,kbuff, fluxl[iter-myBf.begin()]);
                    iter->BufferRecv[4+lr].setValue(i,j,kbuff, fluxr[iter-myBf.begin()]);
                }
            }
            
        }
    }
}


void nuc3d::boundaryCondition::setVisBC(PDEData3d &myPDE,
                                        physicsModel &myPhyMod,
                                        NaiverStokesData3d &myFluxes,
                                        VectorBuffer &myBf)
{
    
    for(auto iter=BCTopo.begin();iter!=BCTopo.end();iter++)
    {
        int iface=static_cast<int>(iter-BCTopo.begin());
        if (((-1)==iter->Type))
        {
            (this->*myVisSetter[BCTopo[iface].id])(myPDE,myPhyMod,myFluxes,myBf,iface);
        }
    }
    
}

void nuc3d::boundaryCondition::setVisBC_Inlet(PDEData3d &myPDE,
                                              physicsModel &myPhyMod,
                                              NaiverStokesData3d &myFluxes,
                                              VectorBuffer &myBf,
                                              int iface)
{
    switch (iface) {
        case 0:
            VisBCsetter_inlet_xi(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 1:
            VisBCsetter_inlet_xi(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 2:
            VisBCsetter_inlet_eta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 3:
            VisBCsetter_inlet_eta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 4:
            VisBCsetter_inlet_zeta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 5:
            VisBCsetter_inlet_zeta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
    }
    
}

void nuc3d::boundaryCondition::VisBCsetter_inlet_xi(PDEData3d &myPDE,
                                                    physicsModel &myPhyMod,
                                                    NaiverStokesData3d &myFluxes,
                                                    VectorBuffer &myBf,
                                                    int lr)
{
    std::vector<double> &q0=BCvalue[lr];
    
    double rho=q0[0];
    double u=q0[1];
    double v=q0[2];
    double w=q0[3];
    double p=q0[4];
    double T;
    double e;
    double alpha;
    
    myPhyMod.prim2conPoint(rho,u,v,w,p,T,e,alpha);
    
    std::vector<double> value={u,v,w,T};
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(auto iter=value.begin();iter!=value.end();iter++)
    {
        for(int k=0;k<nz;k++)
        {
            for(int j=0;j<ny;j++)
            {
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    myBf[iter-value.begin()].BufferRecv[lr].setValue(ibf, j, k,*iter);
                }
            }
        }
    }
}


void nuc3d::boundaryCondition::VisBCsetter_inlet_eta(PDEData3d &myPDE,
                                                     physicsModel &myPhyMod,
                                                     NaiverStokesData3d &myFluxes,
                                                     VectorBuffer &myBf,
                                                     int lr)
{
    std::vector<double> &q0=BCvalue[2+lr];
    
    double rho=q0[0];
    double u=q0[1];
    double v=q0[2];
    double w=q0[3];
    double p=q0[4];
    double T;
    double e;
    double alpha;
    
    myPhyMod.prim2conPoint(rho,u,v,w,p,T,e,alpha);
    
    std::vector<double> value={u,v,w,T};
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(auto iter=value.begin();iter!=value.end();iter++)
    {
        
        for(int k=0;k<nz;k++)
        {
            for(int i=0;i<nx;i++)
            {
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    myBf[iter-value.begin()].BufferRecv[2+lr].setValue(i,ibf, k, *iter);
                }
            }
        }
    }
}

void nuc3d::boundaryCondition::VisBCsetter_inlet_zeta(PDEData3d &myPDE,
                                                      physicsModel &myPhyMod,
                                                      NaiverStokesData3d &myFluxes,
                                                      VectorBuffer &myBf,
                                                      int lr)
{
    std::vector<double> &q0=BCvalue[4+lr];
    
    double rho=q0[0];
    double u=q0[1];
    double v=q0[2];
    double w=q0[3];
    double p=q0[4];
    double T;
    double e;
    double alpha;
    
    myPhyMod.prim2conPoint(rho,u,v,w,p,T,e,alpha);
    
    std::vector<double> value={u,v,w,T};
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(auto iter=value.begin();iter!=value.end();iter++)
    {
        
        for(int j=0;j<ny;j++)
        {
            for(int i=0;i<nx;i++)
            {
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    myBf[iter-value.begin()].BufferRecv[4+lr].setValue(i,j,ibf, *iter);
                }
            }
        }
    }
}

// Outlet condition
void nuc3d::boundaryCondition::setVisBC_Outlet(PDEData3d &myPDE,
                                               physicsModel &myPhyMod,
                                               NaiverStokesData3d &myFluxes,
                                               VectorBuffer &myBf,
                                               int iface)
{
    switch (iface) {
        case 0:
            VisBCsetter_outlet_xi(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 1:
            VisBCsetter_outlet_xi(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 2:
            VisBCsetter_outlet_eta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 3:
            VisBCsetter_outlet_eta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 4:
            VisBCsetter_outlet_zeta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 5:
            VisBCsetter_outlet_zeta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
    }
}

void nuc3d::boundaryCondition::VisBCsetter_outlet_xi(PDEData3d &myPDE,
                                                     physicsModel &myPhyMod,
                                                     NaiverStokesData3d &myFluxes,
                                                     VectorBuffer &myBf,
                                                     int lr)
{
    Field &u=(myFluxes.getPrimatives())[1];
    Field &v=(myFluxes.getPrimatives())[2];
    Field &w=(myFluxes.getPrimatives())[3];
    Field &T=(myFluxes.getAcoustics())[0];
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(int k=0;k<nz;k++)
    {
        for(int j=0;j<ny;j++)
        {
            double data=u.getValue(lr*(nx-1),j,k);
            
            for(int ibf=0;ibf<bfsize;ibf++)
            {
                myBf[0].BufferRecv[lr].setValue(ibf, j, k, data);
            }
        }
    }
    
    for(int k=0;k<nz;k++)
    {
        for(int j=0;j<ny;j++)
        {
            double data=v.getValue(lr*(nx-1),j,k);
            
            for(int ibf=0;ibf<bfsize;ibf++)
            {
                myBf[1].BufferRecv[lr].setValue(ibf, j, k, data);
            }
        }
    }
    
    for(int k=0;k<nz;k++)
    {
        for(int j=0;j<ny;j++)
        {
            double data=w.getValue(lr*(nx-1),j,k);
            
            for(int ibf=0;ibf<bfsize;ibf++)
            {
                myBf[2].BufferRecv[lr].setValue(ibf, j, k, data);
            }
        }
    }
    
    for(int k=0;k<nz;k++)
    {
        for(int j=0;j<ny;j++)
        {
            double data=T.getValue(lr*(nx-1),j,k);
            
            for(int ibf=0;ibf<bfsize;ibf++)
            {
                myBf[3].BufferRecv[lr].setValue(ibf, j, k, data);
            }
        }
    }
}

void nuc3d::boundaryCondition::VisBCsetter_outlet_eta(PDEData3d &myPDE,
                                                      physicsModel &myPhyMod,
                                                      NaiverStokesData3d &myFluxes,
                                                      VectorBuffer &myBf,
                                                      int lr)
{
    Field &u=(myFluxes.getPrimatives())[1];
    Field &v=(myFluxes.getPrimatives())[2];
    Field &w=(myFluxes.getPrimatives())[3];
    Field &T=(myFluxes.getAcoustics())[0];
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(int k=0;k<nz;k++)
    {
        for(int i=0;i<nx;i++)
        {
            double data=u.getValue(i,lr*(ny-1),k);
            
            for(int ibf=0;ibf<bfsize;ibf++)
            {
                myBf[0].BufferRecv[2+lr].setValue(i,ibf, k, data);
            }
        }
    }
    
    for(int k=0;k<nz;k++)
    {
        for(int i=0;i<nx;i++)
        {
            double data=v.getValue(i,lr*(ny-1),k);
            
            for(int ibf=0;ibf<bfsize;ibf++)
            {
                myBf[1].BufferRecv[2+lr].setValue(i,ibf, k, data);
            }
        }
    }
    
    for(int k=0;k<nz;k++)
    {
        for(int i=0;i<nx;i++)
        {
            double data=w.getValue(i,lr*(ny-1),k);
            
            for(int ibf=0;ibf<bfsize;ibf++)
            {
                myBf[2].BufferRecv[2+lr].setValue(i,ibf, k, data);
            }
        }
    }
    
    for(int k=0;k<nz;k++)
    {
        for(int i=0;i<nx;i++)
        {
            double data=T.getValue(i,lr*(ny-1),k);
            
            for(int ibf=0;ibf<bfsize;ibf++)
            {
                myBf[3].BufferRecv[2+lr].setValue(i,ibf, k, data);
            }
        }
    }
}

void nuc3d::boundaryCondition::VisBCsetter_outlet_zeta(PDEData3d &myPDE,
                                                       physicsModel &myPhyMod,
                                                       NaiverStokesData3d &myFluxes,
                                                       VectorBuffer &myBf,
                                                       int lr)
{
    Field &u=(myFluxes.getPrimatives())[1];
    Field &v=(myFluxes.getPrimatives())[2];
    Field &w=(myFluxes.getPrimatives())[3];
    Field &T=(myFluxes.getAcoustics())[0];
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(int j=0;j<ny;j++)
    {
        for(int i=0;i<nx;i++)
        {
            double data=u.getValue(i,j,lr*(nz-1));
            for(int ibf=0;ibf<bfsize;ibf++)
            {
                myBf[0].BufferRecv[4+lr].setValue(i,j,ibf, data);
            }
        }
    }
    
    for(int j=0;j<ny;j++)
    {
        for(int i=0;i<nx;i++)
        {
            double data=v.getValue(i,j,lr*(nz-1));
            for(int ibf=0;ibf<bfsize;ibf++)
            {
                myBf[1].BufferRecv[4+lr].setValue(i,j,ibf, data);
            }
        }
    }
    
    for(int j=0;j<ny;j++)
    {
        for(int i=0;i<nx;i++)
        {
            double data=w.getValue(i,j,lr*(nz-1));
            for(int ibf=0;ibf<bfsize;ibf++)
            {
                myBf[2].BufferRecv[4+lr].setValue(i,j,ibf, data);
            }
        }
    }
    
    for(int j=0;j<ny;j++)
    {
        for(int i=0;i<nx;i++)
        {
            double data=T.getValue(i,j,lr*(nz-1));
            for(int ibf=0;ibf<bfsize;ibf++)
            {
                myBf[3].BufferRecv[4+lr].setValue(i,j,ibf, data);
            }
        }
    }
}

void nuc3d::boundaryCondition::setVisBC_wall(PDEData3d &myPDE,
                                             physicsModel &myPhyMod,
                                             NaiverStokesData3d &myFluxes,
                                             VectorBuffer &myBf,
                                             int iface)
{
    switch (iface) {
        case 0:
            VisBCsetter_wall_xi(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 1:
            VisBCsetter_wall_xi(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 2:
            VisBCsetter_wall_eta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 3:
            VisBCsetter_wall_eta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 4:
            VisBCsetter_wall_zeta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 5:
            VisBCsetter_wall_zeta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
    }
}

void nuc3d::boundaryCondition::VisBCsetter_wall_xi(PDEData3d &myPDE,
                                                   physicsModel &myPhyMod,
                                                   NaiverStokesData3d &myFluxes,
                                                   VectorBuffer &myBf,
                                                   int lr)
{
    std::vector<double> q(4);
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    double T_wall=myPhyMod.getWallTemp();
    
    q[0]=0.0;
    q[1]=0.0;
    q[2]=0.0;
    q[3]=T_wall;
    
    for(auto iter=q.begin();iter!=q.end();iter++)
    {
        for(int k=0;k<nz;k++)
        {
            for(int j=0;j<ny;j++)
            {
                for(int i=0;i<bfsize;i++)
                {
                    myBf[iter-q.begin()].BufferRecv[lr].setValue(i, j, k, *iter);
                }
            }
            
        }
    }
}

void nuc3d::boundaryCondition::VisBCsetter_wall_eta(PDEData3d &myPDE,
                                                    physicsModel &myPhyMod,
                                                    NaiverStokesData3d &myFluxes,
                                                    VectorBuffer &myBf,
                                                    int lr)
{
    std::vector<double> q(4);
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    double T_wall=myPhyMod.getWallTemp();
    
    q[0]=0.0;
    q[1]=0.0;
    q[2]=0.0;
    q[3]=T_wall;
    
    for(auto iter=q.begin();iter!=q.end();iter++)
    {
        for(int k=0;k<nz;k++)
        {
            for(int j=0;j<bfsize;j++)
            {
                for(int i=0;i<nx;i++)
                {
                    myBf[iter-q.begin()].BufferRecv[2+lr].setValue(i,j, k, *iter);
                }
            }
            
        }
    }
}

void nuc3d::boundaryCondition::VisBCsetter_wall_zeta(PDEData3d &myPDE,
                                                     physicsModel &myPhyMod,
                                                     NaiverStokesData3d &myFluxes,
                                                     VectorBuffer &myBf,
                                                     int lr)
{
    std::vector<double> q(4);
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    double T_wall=myPhyMod.getWallTemp();
    
    q[0]=0.0;
    q[1]=0.0;
    q[2]=0.0;
    q[3]=T_wall;
    
    
    for(auto iter=q.begin();iter!=q.end();iter++)
    {        
        for(int k=0;k<bfsize;k++)
        {
            for(int j=0;j<ny;j++)
            {
                for(int i=0;i<nx;i++)
                {
                    myBf[iter-q.begin()].BufferRecv[4+lr].setValue(i,j,k, *iter);
                }
            }
            
        }
    }
}


void nuc3d::boundaryCondition::setVisBC_symm(PDEData3d &myPDE,
                                             physicsModel &myPhyMod,
                                             NaiverStokesData3d &myFluxes,
                                             VectorBuffer &myBf,
                                             int iface)
{
    switch (iface) {
        case 0:
            VisBCsetter_symm_xi(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 1:
            VisBCsetter_symm_xi(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 2:
            VisBCsetter_symm_eta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 3:
            VisBCsetter_symm_eta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 4:
            VisBCsetter_symm_zeta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 5:
            VisBCsetter_symm_zeta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
    }
    
}

void nuc3d::boundaryCondition::VisBCsetter_symm_xi(PDEData3d &myPDE,
                                                   physicsModel &myPhyMod,
                                                   NaiverStokesData3d &myFluxes,
                                                   VectorBuffer &myBf,
                                                   int lr)
{
    std::vector<double> q(4);
    
    VectorField &prim=myFluxes.getPrimatives();
    VectorField &accu=myFluxes.getAcoustics();
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(int k=0;k<nz;k++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int i=0;i<bfsize;i++)
            {
                int iblock=lr*(nx-1)+(1-2*lr)*i;
                int ibuff=(bfsize-1)*(1-lr)-(1-2*lr)*i;
                
                q[0]=prim[1].getValue(iblock,j,k);
                q[1]=prim[2].getValue(iblock,j,k);
                q[2]=prim[3].getValue(iblock,j,k);
                q[3]=accu[0].getValue(iblock,j,k);
                
                for(auto iter=q.begin();iter!=q.end();iter++)
                {
                    myBf[iter-q.begin()].BufferRecv[lr].setValue(ibuff, j, k, *iter);
                }
            }
            
        }
    }
}

void nuc3d::boundaryCondition::VisBCsetter_symm_eta(PDEData3d &myPDE,
                                                    physicsModel &myPhyMod,
                                                    NaiverStokesData3d &myFluxes,
                                                    VectorBuffer &myBf,
                                                    int lr)
{
    std::vector<double> q(4);
    
    VectorField &prim=myFluxes.getPrimatives();
    VectorField &accu=myFluxes.getAcoustics();
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(int k=0;k<nz;k++)
    {
        for(int j=0;j<bfsize;j++)
        {
            for(int i=0;i<nx;i++)
            {
                int jblock=lr*(ny-1)+(1-2*lr)*j;
                int jbuff=(bfsize-1)*(1-lr)-(1-2*lr)*j;
                
                q[0]=prim[1].getValue(i,jblock,k);
                q[1]=prim[2].getValue(i,jblock,k);
                q[2]=prim[3].getValue(i,jblock,k);
                q[3]=accu[0].getValue(i,jblock,k);
                
                for(auto iter=q.begin();iter!=q.end();iter++)
                {
                    myBf[iter-q.begin()].BufferRecv[2+lr].setValue(i,jbuff, k, *iter);
                }
            }
            
        }
    }
}

void nuc3d::boundaryCondition::VisBCsetter_symm_zeta(PDEData3d &myPDE,
                                                     physicsModel &myPhyMod,
                                                     NaiverStokesData3d &myFluxes,
                                                     VectorBuffer &myBf,
                                                     int lr)
{
    std::vector<double> q(4);
    
    VectorField &prim=myFluxes.getPrimatives();
    VectorField &accu=myFluxes.getAcoustics();
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    for(int k=0;k<bfsize;k++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int i=0;i<nx;i++)
            {
                int kblock=lr*(nz-1)+(1-2*lr)*k;
                int kbuff=(bfsize-1)*(1-lr)-(1-2*lr)*k;
                q[0]=prim[1].getValue(i,j,kblock);
                q[1]=prim[2].getValue(i,j,kblock);
                q[2]=prim[3].getValue(i,j,kblock);
                q[3]=accu[0].getValue(i,j,kblock);
                
                for(auto iter=q.begin();iter!=q.end();iter++)
                {
                    myBf[iter-q.begin()].BufferRecv[4+lr].setValue(i,j,kbuff, *iter);
                }
            }
            
        }
    }
}

void nuc3d::boundaryCondition::setVisFluxBC(PDEData3d &myPDE,
                                            physicsModel &myPhyMod,
                                            NaiverStokesData3d &myFluxes,
                                            VectorBuffer &myBf)
{
    
    for(auto iter=BCTopo.begin();iter!=BCTopo.end();iter++)
    {
        int iface=static_cast<int>(iter-BCTopo.begin());
        if (((-1)==iter->Type))
        {
            (this->*myVisFluxSetter[BCTopo[iface].id])(myPDE,myPhyMod,myFluxes,myBf,iface);
        }
    }
    
}

void nuc3d::boundaryCondition::setVisFluxBC_Inlet(PDEData3d &myPDE,
                                                  physicsModel &myPhyMod,
                                                  NaiverStokesData3d &myFluxes,
                                                  VectorBuffer &myBf,
                                                  int iface)
{
    switch (iface) {
        case 0:
            VisFluxBCsetter_inlet_xi(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 1:
            VisFluxBCsetter_inlet_xi(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 2:
            VisFluxBCsetter_inlet_eta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 3:
            VisFluxBCsetter_inlet_eta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 4:
            VisFluxBCsetter_inlet_zeta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 5:
            VisFluxBCsetter_inlet_zeta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
    }
    
}

void nuc3d::boundaryCondition::VisFluxBCsetter_inlet_xi(PDEData3d &myPDE,
                                                        physicsModel &myPhyMod,
                                                        NaiverStokesData3d &myFluxes,
                                                        VectorBuffer &myBf,
                                                        int lr)
{
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    VectorField &visflux_xi=myFluxes.getVisFlux_xi();
    
    for(auto iter=visflux_xi.begin();iter!=visflux_xi.end();iter++)
    {
        for(int k=0;k<nz;k++)
        {
            for(int j=0;j<ny;j++)
            {
                for(int i=0;i<bfsize;i++)
                {
                    double value=0.0;
                    
                    myBf[iter-visflux_xi.begin()].BufferRecv[lr].setValue(i, j, k, value);
                }
            }
            
        }
    }
}

void nuc3d::boundaryCondition::VisFluxBCsetter_inlet_eta(PDEData3d &myPDE,
                                                         physicsModel &myPhyMod,
                                                         NaiverStokesData3d &myFluxes,
                                                         VectorBuffer &myBf,
                                                         int lr)
{
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    VectorField &visflux_eta=myFluxes.getVisFlux_eta();
    
    for(auto iter=visflux_eta.begin();iter!=visflux_eta.end();iter++)
    {
        
        for(int k=0;k<nz;k++)
        {
            for(int j=0;j<bfsize;j++)
            {
                for(int i=0;i<nx;i++)
                {
                    double value=0.0;
                    myBf[iter-visflux_eta.begin()].BufferRecv[2+lr].setValue(i,j, k, value);
                }
            }
            
        }
    }
    
}

void nuc3d::boundaryCondition::VisFluxBCsetter_inlet_zeta(PDEData3d &myPDE,
                                                          physicsModel &myPhyMod,
                                                          NaiverStokesData3d &myFluxes,
                                                          VectorBuffer &myBf,
                                                          int lr)
{
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    VectorField &visflux_zeta=myFluxes.getVisFlux_zeta();
    
    for(auto iter=visflux_zeta.begin();iter!=visflux_zeta.end();iter++)
    {
        for(int k=0;k<bfsize;k++)
        {
            for(int j=0;j<ny;j++)
            {
                for(int i=0;i<nx;i++)
                {
                    double value=0.0;
                    
                    myBf[iter-visflux_zeta.begin()].BufferRecv[4+lr].setValue(i,j,k, value);
                }
            }
        }
    }
    
}


//outlet
void nuc3d::boundaryCondition::setVisFluxBC_Outlet(PDEData3d &myPDE,
                                                   physicsModel &myPhyMod,
                                                   NaiverStokesData3d &myFluxes,
                                                   VectorBuffer &myBf,
                                                   int iface)
{
    switch (iface) {
        case 0:
            VisFluxBCsetter_outlet_xi(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 1:
            VisFluxBCsetter_outlet_xi(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 2:
            VisFluxBCsetter_outlet_eta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 3:
            VisFluxBCsetter_outlet_eta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 4:
            VisFluxBCsetter_outlet_zeta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 5:
            VisFluxBCsetter_outlet_zeta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
    }
    
    
}

void nuc3d::boundaryCondition::VisFluxBCsetter_outlet_xi(PDEData3d &myPDE,
                                                         physicsModel &myPhyMod,
                                                         NaiverStokesData3d &myFluxes,
                                                         VectorBuffer &myBf,
                                                         int lr)
{
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    VectorField &visflux_xi=myFluxes.getVisFlux_xi();
    
    for(auto iter=visflux_xi.begin();iter!=visflux_xi.end();iter++)
    {
        for(int k=0;k<nz;k++)
        {
            for(int j=0;j<ny;j++)
            {
                for(int i=0;i<bfsize;i++)
                {
                    double value=iter->getValue(lr*(nx-1), j, k);
                    
                    myBf[iter-visflux_xi.begin()].BufferRecv[lr].setValue(i, j, k, value);
                }
            }
            
        }
    }
    
    
}
void nuc3d::boundaryCondition::VisFluxBCsetter_outlet_eta(PDEData3d &myPDE,
                                                          physicsModel &myPhyMod,
                                                          NaiverStokesData3d &myFluxes,
                                                          VectorBuffer &myBf,
                                                          int lr)
{
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    VectorField &visflux_eta=myFluxes.getVisFlux_eta();
    
    for(auto iter=visflux_eta.begin();iter!=visflux_eta.end();iter++)
    {
        
        for(int k=0;k<nz;k++)
        {
            for(int j=0;j<bfsize;j++)
            {
                for(int i=0;i<nx;i++)
                {
                    double value=iter->getValue(i, lr*(ny-1), k);
                    myBf[iter-visflux_eta.begin()].BufferRecv[2+lr].setValue(i,j, k, value);
                }
            }
            
        }
    }
    
    
}

void nuc3d::boundaryCondition::VisFluxBCsetter_outlet_zeta(PDEData3d &myPDE,
                                                           physicsModel &myPhyMod,
                                                           NaiverStokesData3d &myFluxes,
                                                           VectorBuffer &myBf,
                                                           int lr)
{
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    VectorField &visflux_zeta=myFluxes.getVisFlux_zeta();
    
    for(auto iter=visflux_zeta.begin();iter!=visflux_zeta.end();iter++)
    {
        for(int k=0;k<bfsize;k++)
        {
            for(int j=0;j<ny;j++)
            {
                for(int i=0;i<nx;i++)
                {
                    double value=iter->getValue(i, j, lr*(nz-1));
                    
                    myBf[iter-visflux_zeta.begin()].BufferRecv[4+lr].setValue(i,j,k, value);
                }
            }
        }
    }
    
    
}

//wall
void nuc3d::boundaryCondition::setVisFluxBC_wall(PDEData3d &myPDE,
                                                 physicsModel &myPhyMod,
                                                 NaiverStokesData3d &myFluxes,
                                                 VectorBuffer &myBf,
                                                 int iface)
{
    switch (iface) {
        case 0:
            VisFluxBCsetter_wall_xi(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 1:
            VisFluxBCsetter_wall_xi(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 2:
            VisFluxBCsetter_wall_eta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 3:
            VisFluxBCsetter_wall_eta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 4:
            VisFluxBCsetter_wall_zeta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 5:
            VisFluxBCsetter_wall_zeta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
    }
    
    
}

void nuc3d::boundaryCondition::VisFluxBCsetter_wall_xi(PDEData3d &myPDE,
                                                       physicsModel &myPhyMod,
                                                       NaiverStokesData3d &myFluxes,
                                                       VectorBuffer &myBf,
                                                       int lr)
{
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    VectorField &visflux_xi=myFluxes.getVisFlux_xi();
    
    for(auto iter=visflux_xi.begin();iter!=visflux_xi.end();iter++)
    {
        for(int k=0;k<nz;k++)
        {
            for(int j=0;j<ny;j++)
            {
                for(int i=0;i<bfsize;i++)
                {
                    int iblock=lr*(nx-1)+(1-2*lr)*i;
                    int ibuff=(bfsize-1)*(1-lr)-(1-2*lr)*i;
                    
                    double value=iter->getValue(iblock, j, k);
                    
                    myBf[iter-visflux_xi.begin()].BufferRecv[lr].setValue(ibuff, j, k, value);
                }
            }
            
        }
    }
    
    
}

void nuc3d::boundaryCondition::VisFluxBCsetter_wall_eta(PDEData3d &myPDE,
                                                        physicsModel &myPhyMod,
                                                        NaiverStokesData3d &myFluxes,
                                                        VectorBuffer &myBf,
                                                        int lr)
{
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    VectorField &visflux_eta=myFluxes.getVisFlux_eta();
    
    for(auto iter=visflux_eta.begin();iter!=visflux_eta.end();iter++)
    {
        
        for(int k=0;k<nz;k++)
        {
            for(int j=0;j<bfsize;j++)
            {
                for(int i=0;i<nx;i++)
                {
                    int jblock=lr*(ny-1)+(1-2*lr)*j;
                    int jbuff=(bfsize-1)*(1-lr)-(1-2*lr)*j;
                    
                    double value=iter->getValue(i, jblock, k);
                    myBf[iter-visflux_eta.begin()].BufferRecv[2+lr].setValue(i,jbuff, k, value);
                }
            }
            
        }
    }
    
    
}
void nuc3d::boundaryCondition::VisFluxBCsetter_wall_zeta(PDEData3d &myPDE,
                                                         physicsModel &myPhyMod,
                                                         NaiverStokesData3d &myFluxes,
                                                         VectorBuffer &myBf,
                                                         int lr)
{
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    VectorField &visflux_zeta=myFluxes.getVisFlux_zeta();
    
    for(auto iter=visflux_zeta.begin();iter!=visflux_zeta.end();iter++)
    {
        for(int k=0;k<bfsize;k++)
        {
            for(int j=0;j<ny;j++)
            {
                for(int i=0;i<nx;i++)
                {
                    int kblock=lr*(nz-1)+(1-2*lr)*k;
                    int kbuff=(bfsize-1)*(1-lr)-(1-2*lr)*k;
                    
                    double value=iter->getValue(i, j, kblock);
                    
                    
                    myBf[iter-visflux_zeta.begin()].BufferRecv[4+lr].setValue(i,j,kbuff, value);
                }
            }
        }
    }
    
    
}

//symmetric
void nuc3d::boundaryCondition::setVisFluxBC_symm(PDEData3d &myPDE,
                                                 physicsModel &myPhyMod,
                                                 NaiverStokesData3d &myFluxes,
                                                 VectorBuffer &myBf,
                                                 int iface)
{
    
    switch (iface) {
        case 0:
            VisFluxBCsetter_symm_xi(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 1:
            VisFluxBCsetter_symm_xi(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 2:
            VisFluxBCsetter_symm_eta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 3:
            VisFluxBCsetter_symm_eta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
        case 4:
            VisFluxBCsetter_symm_zeta(myPDE,myPhyMod,myFluxes,myBf,0);
            break;
        case 5:
            VisFluxBCsetter_symm_zeta(myPDE,myPhyMod,myFluxes,myBf,1);
            break;
    }
    
}

void nuc3d::boundaryCondition::VisFluxBCsetter_symm_xi(PDEData3d &myPDE,
                                                       physicsModel &myPhyMod,
                                                       NaiverStokesData3d &myFluxes,
                                                       VectorBuffer &myBf,
                                                       int lr)
{
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    VectorField &visflux_xi=myFluxes.getVisFlux_xi();
    
    for(auto iter=visflux_xi.begin();iter!=visflux_xi.end();iter++)
    {
        for(int k=0;k<nz;k++)
        {
            for(int j=0;j<ny;j++)
            {
                for(int i=0;i<bfsize;i++)
                {
                    int iblock=lr*(nx-1)+(1-2*lr)*i;
                    int ibuff=(bfsize-1)*(1-lr)-(1-2*lr)*i;
                    
                    double value=iter->getValue(iblock, j, k);
                    
                    myBf[iter-visflux_xi.begin()].BufferRecv[lr].setValue(ibuff, j, k, value);
                }
            }
            
        }
    }
    
    
}
void nuc3d::boundaryCondition::VisFluxBCsetter_symm_eta(PDEData3d &myPDE,
                                                        physicsModel &myPhyMod,
                                                        NaiverStokesData3d &myFluxes,
                                                        VectorBuffer &myBf,
                                                        int lr)
{
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    VectorField &visflux_eta=myFluxes.getVisFlux_eta();
    
    for(auto iter=visflux_eta.begin();iter!=visflux_eta.end();iter++)
    {
        
        for(int k=0;k<nz;k++)
        {
            for(int j=0;j<bfsize;j++)
            {
                for(int i=0;i<nx;i++)
                {
                    int jblock=lr*(ny-1)+(1-2*lr)*j;
                    int jbuff=(bfsize-1)*(1-lr)-(1-2*lr)*j;
                    
                    double value=iter->getValue(i, jblock, k);
                    myBf[iter-visflux_eta.begin()].BufferRecv[2+lr].setValue(i,jbuff, k, value);
                }
            }
            
        }
    }
    
    
}
void nuc3d::boundaryCondition::VisFluxBCsetter_symm_zeta(PDEData3d &myPDE,
                                                         physicsModel &myPhyMod,
                                                         NaiverStokesData3d &myFluxes,
                                                         VectorBuffer &myBf,
                                                         int lr)
{
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    VectorField &visflux_zeta=myFluxes.getVisFlux_zeta();
    
    for(auto iter=visflux_zeta.begin();iter!=visflux_zeta.end();iter++)
    {
        for(int k=0;k<bfsize;k++)
        {
            for(int j=0;j<ny;j++)
            {
                for(int i=0;i<nx;i++)
                {
                    int kblock=lr*(nz-1)+(1-2*lr)*k;
                    int kbuff=(bfsize-1)*(1-lr)-(1-2*lr)*k;
                    
                    double value=iter->getValue(i, j, kblock);
                    
                    
                    myBf[iter-visflux_zeta.begin()].BufferRecv[4+lr].setValue(i,j,kbuff, value);
                }
            }
        }
    }
}


