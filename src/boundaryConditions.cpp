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
                double *pSend=iter->BufferSend[lr].getDataPtr();
                double *pRecv=iter->BufferRecv[lr].getDataPtr();
                
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    int idx_xi=bfsize*ny*k+bfsize*j+ibf;
                    
                    pSend[idx_xi]=fluxl[iter-myBf.begin()];
                    pRecv[idx_xi]=fluxr[iter-myBf.begin()];
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
                double *pSend=iter->BufferSend[2+lr].getDataPtr();
                double *pRecv=iter->BufferRecv[2+lr].getDataPtr();
                
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    int idx_xi=bfsize*nz*i+bfsize*k+ibf;
                    
                    pSend[idx_xi]=fluxl[iter-myBf.begin()];
                    pRecv[idx_xi]=fluxr[iter-myBf.begin()];
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
                double *pSend=iter->BufferSend[4+lr].getDataPtr();
                double *pRecv=iter->BufferRecv[4+lr].getDataPtr();
                
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    int idx_xi=bfsize*nx*j+bfsize*i+ibf;
                    
                    pSend[idx_xi]=fluxl[iter-myBf.begin()];
                    pRecv[idx_xi]=fluxr[iter-myBf.begin()];
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
                double *pSend=iter->BufferSend[lr].getDataPtr();
                double *pRecv=iter->BufferRecv[lr].getDataPtr();
                
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    int idx_xi=bfsize*ny*k+bfsize*j+ibf;
                    
                    pSend[idx_xi]=fluxl[iter-myBf.begin()];
                    pRecv[idx_xi]=fluxr[iter-myBf.begin()];
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
                double *pSend=iter->BufferSend[2+lr].getDataPtr();
                double *pRecv=iter->BufferRecv[2+lr].getDataPtr();
                
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    int idx_xi=bfsize*nz*i+bfsize*k+ibf;
                    
                    pSend[idx_xi]=fluxl[iter-myBf.begin()];
                    pRecv[idx_xi]=fluxr[iter-myBf.begin()];
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
    
    for(int j=0;j<ny;j++)
    {
        for(int i=0;i<nx;i++)
        {
            double jacob=jac.getValue(i,j,lr*(nz-1));
            double zeta_x=zeta_xyz[0].getValue(i,j,lr*(nz-1));
            double zeta_y=zeta_xyz[1].getValue(i,j,lr*(nz-1));
            double zeta_z=zeta_xyz[2].getValue(i,j,lr*(nz-1));
            
            for(auto iter=q.begin();iter!=q.end();iter++)
                *iter=prim[iter-q.begin()].getValue(i,j,lr*(nz-1));
            
            myPhyMod.solveRiemannPoint(q, jacob, zeta_x, zeta_y, zeta_z, fluxl, fluxr);
            
            for(auto iter=myBf.begin();iter!=myBf.end();iter++)
            {
                double *pSend=iter->BufferSend[4+lr].getDataPtr();
                double *pRecv=iter->BufferRecv[4+lr].getDataPtr();
                
                for(int ibf=0;ibf<bfsize;ibf++)
                {
                    int idx_xi=bfsize*nx*j+bfsize*i+ibf;
                    
                    pSend[idx_xi]=fluxl[iter-myBf.begin()];
                    pRecv[idx_xi]=fluxr[iter-myBf.begin()];
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
    VectorField &eta_xyz=myFluxes.getEta_xyz();
    VectorField &zeta_xyz=myFluxes.getZeta_xyz();
    
    VectorField &dx=myFluxes.getDx();
    VectorField &dy=myFluxes.getDy();
    VectorField &dz=myFluxes.getDz();
    
    VectorField &prim=myFluxes.getPrimatives();
    
    double t_wall=myPhyMod.getWallTemp();
    double gamma=myPhyMod.getGamma();
    double mach=myPhyMod.getMach();
    
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
                //
                //                q[1]=-x_xi*U-x_eta*V-x_zeta*W;
                //                q[2]=-y_xi*U-y_eta*V-y_zeta*W;
                //                q[3]=-z_xi*U-z_eta*V-z_zeta*W;
                
                q[0]=q[4]*gamma*mach*mach/t_wall;
                q[1]=0.0;
                q[2]=0.0;
                q[3]=0.0;
                
                
                myPhyMod.solveRiemannPoint(q, jacob, xi_x, xi_y, xi_z, fluxl, fluxr);
                
                int idx_xi=bfsize*ny*k+bfsize*j+ibuff;
                
                for(auto iter=myBf.begin();iter!=myBf.end();iter++)
                {
                    double *pSend=iter->BufferSend[lr].getDataPtr();
                    double *pRecv=iter->BufferRecv[lr].getDataPtr();
                    
                    pSend[idx_xi]=fluxl[iter-myBf.begin()];
                    pRecv[idx_xi]=fluxr[iter-myBf.begin()];
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
    VectorField &xi_xyz=myFluxes.getXi_xyz();
    VectorField &eta_xyz=myFluxes.getEta_xyz();
    VectorField &zeta_xyz=myFluxes.getZeta_xyz();
    
    VectorField &dx=myFluxes.getDx();
    VectorField &dy=myFluxes.getDy();
    VectorField &dz=myFluxes.getDz();
    
    VectorField &prim=myFluxes.getPrimatives();
    
    double t_wall=myPhyMod.getWallTemp();
    double gamma=myPhyMod.getGamma();
    double mach=myPhyMod.getMach();

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
                //
                //                q[1]=-x_xi*U-x_eta*V-x_zeta*W;
                //                q[2]=-y_xi*U-y_eta*V-y_zeta*W;
                //                q[3]=-z_xi*U-z_eta*V-z_zeta*W;
                
                q[0]=q[4]*gamma*mach*mach/t_wall;
                q[1]=0.0;
                q[2]=0.0;
                q[3]=0.0;
                
                
                myPhyMod.solveRiemannPoint(q, jacob, eta_x, eta_y, eta_z, fluxl, fluxr);
                
                int idx_xi=bfsize*nz*i+bfsize*k+jbuff;
                
                for(auto iter=myBf.begin();iter!=myBf.end();iter++)
                {
                    double *pSend=iter->BufferSend[2+lr].getDataPtr();
                    double *pRecv=iter->BufferRecv[2+lr].getDataPtr();
                    
                    pSend[idx_xi]=fluxl[iter-myBf.begin()];
                    pRecv[idx_xi]=fluxr[iter-myBf.begin()];
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
    VectorField &xi_xyz=myFluxes.getXi_xyz();
    VectorField &eta_xyz=myFluxes.getEta_xyz();
    VectorField &zeta_xyz=myFluxes.getZeta_xyz();
    
    VectorField &dx=myFluxes.getDx();
    VectorField &dy=myFluxes.getDy();
    VectorField &dz=myFluxes.getDz();
    
    VectorField &prim=myFluxes.getPrimatives();
    
    double t_wall=myPhyMod.getWallTemp();
    double gamma=myPhyMod.getGamma();
    double mach=myPhyMod.getMach();

    
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
                //
                //                q[1]=-x_xi*U-x_eta*V-x_zeta*W;
                //                q[2]=-y_xi*U-y_eta*V-y_zeta*W;
                //                q[3]=-z_xi*U-z_eta*V-z_zeta*W;
                
                q[0]=q[4]*gamma*mach*mach/t_wall;
                q[1]=0.0;
                q[2]=0.0;
                q[3]=0.0;
                
                
                myPhyMod.solveRiemannPoint(q, jacob, zeta_x, zeta_y, zeta_z, fluxl, fluxr);
                
                int idx_xi=bfsize*nx*j+bfsize*i+kbuff;
                
                for(auto iter=myBf.begin();iter!=myBf.end();iter++)
                {
                    double *pSend=iter->BufferSend[4+lr].getDataPtr();
                    double *pRecv=iter->BufferRecv[4+lr].getDataPtr();
                    
                    pSend[idx_xi]=fluxl[iter-myBf.begin()];
                    pRecv[idx_xi]=fluxr[iter-myBf.begin()];
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
                
                int idx_xi=bfsize*ny*k+bfsize*j+ibuff;
                
                for(auto iter=myBf.begin();iter!=myBf.end();iter++)
                {
                    double *pSend=iter->BufferSend[lr].getDataPtr();
                    double *pRecv=iter->BufferRecv[lr].getDataPtr();
                    
                    pSend[idx_xi]=fluxl[iter-myBf.begin()];
                    pRecv[idx_xi]=fluxr[iter-myBf.begin()];
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
                
                int idx_xi=bfsize*nz*i+bfsize*k+jbuff;
                
                for(auto iter=myBf.begin();iter!=myBf.end();iter++)
                {
                    double *pSend=iter->BufferSend[2+lr].getDataPtr();
                    double *pRecv=iter->BufferRecv[2+lr].getDataPtr();
                    
                    pSend[idx_xi]=fluxl[iter-myBf.begin()];
                    pRecv[idx_xi]=fluxr[iter-myBf.begin()];
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
                
                int idx_xi=bfsize*nx*j+bfsize*i+kbuff;
                
                for(auto iter=myBf.begin();iter!=myBf.end();iter++)
                {
                    double *pSend=iter->BufferSend[4+lr].getDataPtr();
                    double *pRecv=iter->BufferRecv[4+lr].getDataPtr();
                    
                    pSend[idx_xi]=fluxl[iter-myBf.begin()];
                    pRecv[idx_xi]=fluxr[iter-myBf.begin()];
                }
            }
            
        }
    }
}
