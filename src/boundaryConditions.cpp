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
    VectorField &flux_xi_l=myFluxes.getFluxXi().FluxL;
    VectorField &flux_xi_r=myFluxes.getFluxXi().FluxR;
    
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    auto begl=fluxl.begin();
    auto endl=fluxl.end();
    
    for(auto iter=begl;iter!=endl;iter++)
    {
        for(int k=0;k<nz;k++)
        {
            for(int j=0;j<ny;j++)
            {
                *iter=flux_xi_l[iter-begl].getValue(lr*(nx-1),j,k);
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    myBf[iter-begl].BufferSend[lr].setValue(ibf, j, k, *iter);
                }
            }
        }
    }
    
    auto begr=fluxr.begin();
    auto endr=fluxr.end();
    
    for(auto iter=begr;iter!=endr;iter++)
    {
        for(int k=0;k<nz;k++)
        {
            for(int j=0;j<ny;j++)
            {
                *iter=flux_xi_r[iter-begr].getValue(lr*(nx-1),j,k);
                
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    myBf[iter-begr].BufferRecv[lr].setValue(ibf, j, k, *iter);
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
    VectorField &flux_eta_l=myFluxes.getFluxEta().FluxL;
    VectorField &flux_eta_r=myFluxes.getFluxEta().FluxR;
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    auto begl=fluxl.begin();
    auto endl=fluxl.end();
    
    for(auto iter=begl;iter!=endl;iter++)
    {
        for(int k=0;k<nz;k++)
        {
            for(int i=0;i<nx;i++)
            {
                *iter=flux_eta_l[iter-begl].getValue(i,lr*(ny-1),k);
                
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    myBf[iter-begl].BufferSend[2+lr].setValue(i,ibf, k, *iter);
                }
            }
        }
    }
    
    auto begr=fluxr.begin();
    auto endr=fluxr.end();
    
    for(auto iter=begr;iter!=endr;iter++)
    {
        for(int k=0;k<nz;k++)
        {
            for(int i=0;i<nx;i++)
            {
                *iter=flux_eta_r[iter-begr].getValue(i,lr*(ny-1),k);
                
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    myBf[iter-begr].BufferRecv[2+lr].setValue(i,ibf, k, *iter);
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
    VectorField &flux_zeta_l=myFluxes.getFluxZeta().FluxL;
    VectorField &flux_zeta_r=myFluxes.getFluxZeta().FluxR;
    
    const int nx=myFluxes.nx;
    const int ny=myFluxes.ny;
    const int nz=myFluxes.nz;
    const int bfsize=myBf[0].bufferWidth;
    
    auto begl=fluxl.begin();
    auto endl=fluxl.end();
    
    for(auto iter=begl;iter!=endl;iter++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int i=0;i<nx;i++)
            {
                *iter=flux_zeta_l[iter-begl].getValue(i,j,lr*(nz-1));
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    myBf[iter-begl].BufferSend[4+lr].setValue(i,j,ibf, *iter);
                }
            }
        }
    }
    
    auto begr=fluxr.begin();
    auto endr=fluxr.end();
    
    for(auto iter=begr;iter!=endr;iter++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int i=0;i<nx;i++)
            {
                *iter=flux_zeta_r[iter-begr].getValue(i,j,lr*(nz-1));
                
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    myBf[iter-begr].BufferRecv[4+lr].setValue(i,j,ibf, *iter);
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
                
                for(auto iter=q.begin();iter!=q.end();iter++)
                    *iter=prim[iter-q.begin()].getValue(iblock, j, k);
                
                //u=-u0 v=-v0 w=-w0
                q[1]*=-1.0;
                q[2]*=-1.0;
                q[3]*=-1.0;
                
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
                double eta_x=eta_xyz[0].getValue(i,jblock,k);
                double eta_y=eta_xyz[1].getValue(i,jblock,k);
                double eta_z=eta_xyz[2].getValue(i,jblock,k);
                
                for(auto iter=q.begin();iter!=q.end();iter++)
                    *iter=prim[iter-q.begin()].getValue(i,jblock, k);
                
                //u=-u0 v=-v0 w=-w0
                q[1]*=-1.0;
                q[2]*=-1.0;
                q[3]*=-1.0;
                
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
                double zeta_x=zeta_xyz[0].getValue(i,j,kblock);
                double zeta_y=zeta_xyz[1].getValue(i,j,kblock);
                double zeta_z=zeta_xyz[2].getValue(i,j,kblock);
                
                for(auto iter=q.begin();iter!=q.end();iter++)
                    *iter=prim[iter-q.begin()].getValue(i,j,kblock);
                
                //u=-u0 v=-v0 w=-w0
                q[1]*=-1.0;
                q[2]*=-1.0;
                q[3]*=-1.0;
                
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
                
                for(auto iter=q.begin();iter!=q.end();iter++)
                    *iter=prim[iter-q.begin()].getValue(iblock, j, k);
                
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
    VectorField &eta_xyz=myFluxes.getEta_xyz();
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
                double eta_x=eta_xyz[0].getValue(i,jblock,k);
                double eta_y=eta_xyz[1].getValue(i,jblock,k);
                double eta_z=eta_xyz[2].getValue(i,jblock,k);
                
                for(auto iter=q.begin();iter!=q.end();iter++)
                    *iter=prim[iter-q.begin()].getValue(i,jblock, k);
                
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
                int kblock=lr*(nz-1)+(1-2*lr)*k;
                int kbuff=(bfsize-1)*(1-lr)-(1-2*lr)*k;
                
                double jacob=jac.getValue(i,j,kblock);
                double zeta_x=zeta_xyz[0].getValue(i,j,kblock);
                double zeta_y=zeta_xyz[1].getValue(i,j,kblock);
                double zeta_z=zeta_xyz[2].getValue(i,j,kblock);
                
                for(auto iter=q.begin();iter!=q.end();iter++)
                    *iter=prim[iter-q.begin()].getValue(i,j,kblock);
                
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
