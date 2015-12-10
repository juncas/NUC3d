//
//  NaiverStokesData3d.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/11/3.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#include "NaiverStokes3d.h"
#include "physicsModel.h"
#include "PDEData3d.hpp"
#include "bufferData.hpp"
nuc3d::gradvector::gradvector(int nx0,int ny0,int nz0):
dxi(nx0,ny0,nz0),
deta(nx0,ny0,nz0),
dzeta(nx0,ny0,nz0)
{
    
}
nuc3d::gradvector::~gradvector()
{}
/*************************************************************************************
 Member functions of class: NaiverStokesData3d
 **************************************************************************************/
nuc3d::NaiverStokesData3d::NaiverStokesData3d(int nx0, int ny0, int nz0, int neqs):
EulerData3D(nx0,ny0,nz0,neqs),
du(nx0,ny0,nz0),
dv(nx0,ny0,nz0),
dw(nx0,ny0,nz0),
dT(nx0,ny0,nz0),
miu(nx0,ny0,nz0,1.0),
coeff(nx0,ny0,nz0,1.0),
tau(9,Field(nx0,ny0,nz0)),
Flux_xi_vis(neqs,Field(nx0,ny0,nz0)),
Flux_eta_vis(neqs,Field(nx0,ny0,nz0)),
Flux_zeta_vis(neqs,Field(nx0,ny0,nz0)),
dfvdxi(neqs,Field(nx0,ny0,nz0)),
dgvdeta(neqs,Field(nx0,ny0,nz0)),
dhvdzeta(neqs,Field(nx0,ny0,nz0))
{
    
}

void nuc3d::NaiverStokesData3d::solve(PDEData3d &myPDE,
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
    
    nuc3d::NaiverStokesData3d::solveVis(myPDE,myOP,myModel,myBf,myMPI,myBC);
    nuc3d::NaiverStokesData3d::solveRHS(myPDE);
}

void nuc3d::NaiverStokesData3d::solveVis(PDEData3d &myPDE,
                                         fieldOperator3d &myOP,
                                         physicsModel &myModel,
                                         std::vector<bufferData> &myBf,
                                         MPIComunicator3d_nonblocking &myMPI,
                                         boundaryCondition &myBC)
{
    setBoundaryGrad(myPDE,myOP,myModel,myBf,myMPI,myBC);
    solveGrads(myPDE, myOP, myBf, myMPI);
    
    solveViscousFlux(myModel);
    setBoundaryViscousFlux();
    
    setDerivativesVis();
}


void nuc3d::NaiverStokesData3d::setBoundaryGrad(PDEData3d &myPDE,
                                                fieldOperator3d &myOP,
                                                physicsModel &myModel,
                                                std::vector<bufferData> &myBf,
                                                MPIComunicator3d_nonblocking &myMPI,
                                                boundaryCondition &myBC)
{
    myBC.setVisBC(myPDE, myModel,*this,myBf);
}

void nuc3d::NaiverStokesData3d::solveGrads(PDEData3d &myPDE,
                                           fieldOperator3d &myOP,
                                           std::vector<bufferData> &myBf,
                                           MPIComunicator3d_nonblocking &myMPI)
{
    Field &u=EulerData3D::W_Euler[1];
    Field &v=EulerData3D::W_Euler[2];
    Field &w=EulerData3D::W_Euler[3];
    Field &T=EulerData3D::W0_Euler[0];
    
    solve_grad(u,myOP,myBf[0],myMPI,du,0);
    solve_grad(v,myOP,myBf[1],myMPI,dv,1);
    solve_grad(w,myOP,myBf[2],myMPI,dw,2);
    solve_grad(T,myOP,myBf[3],myMPI,dT,3);
    MPI_Barrier(MPI_COMM_WORLD);
}

void nuc3d::NaiverStokesData3d::solve_grad(Field &myField,
                                           fieldOperator3d &myOP,
                                           bufferData &myBf,
                                           MPIComunicator3d_nonblocking &myMPI,
                                           gradvector &myGrad,
                                           int fdID)
{
    solveGradXi(myField,myOP,myBf,myMPI,myGrad.getdxi(),fdID);
    solveGradEta(myField,myOP,myBf,myMPI,myGrad.getdeta(),fdID);
    solveGradZeta(myField,myOP,myBf,myMPI,myGrad.getdzeta(),fdID);
}

void nuc3d::NaiverStokesData3d::solveGradXi(Field &myField,
                                            fieldOperator3d &myOP,
                                            bufferData &myBf,
                                            MPIComunicator3d_nonblocking &myMPI,
                                            Field &dxi,
                                            int fdID)
{
    myMPI.bufferSendRecv(myField, myBf, 0, fdID);
    myOP.differenceInner(myField, 0, dxi);
    myMPI.waitSendRecv(myBf, 0);
    myOP.differenceBoundary(myField, myBf.BufferRecv[0], myBf.BufferRecv[1], 0, dxi);
}

void nuc3d::NaiverStokesData3d::solveGradEta(Field &myField,
                                             fieldOperator3d &myOP,
                                             bufferData &myBf,
                                             MPIComunicator3d_nonblocking &myMPI,
                                             Field &deta,
                                             int fdID)
{
    myMPI.bufferSendRecv(myField, myBf, 1, fdID);
    myOP.differenceInner(myField, 1, deta);
    myMPI.waitSendRecv(myBf, 1);
    myOP.differenceBoundary(myField, myBf.BufferRecv[2], myBf.BufferRecv[3], 1,deta);
}

void nuc3d::NaiverStokesData3d::solveGradZeta(Field &myField,
                                              fieldOperator3d &myOP,
                                              bufferData &myBf,
                                              MPIComunicator3d_nonblocking &myMPI,
                                              Field &dzeta,
                                              int fdID)
{
    myMPI.bufferSendRecv(myField, myBf, 2, fdID);
    myOP.differenceInner(myField, 2, dzeta);
    myMPI.waitSendRecv(myBf, 2);
    myOP.differenceBoundary(myField, myBf.BufferRecv[4], myBf.BufferRecv[5], 2, dzeta);
}

void nuc3d::NaiverStokesData3d::solveViscousFlux(physicsModel &myPhyMod)
{
    
    myPhyMod.getMiu(W0_Euler[0], miu, coeff);
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                //u,v,w,T
                double u=W_Euler[1].getValue(i, j, k);
                double v=W_Euler[2].getValue(i, j, k);
                double w=W_Euler[3].getValue(i, j, k);
                
                //Grads of u,v,w,T
                double uxi=du.getdxi().getValue(i, j, k);
                double ueta=du.getdeta().getValue(i, j, k);
                double uzeta=du.getdzeta().getValue(i, j, k);
                
                double vxi=dv.getdxi().getValue(i, j, k);
                double veta=dv.getdeta().getValue(i, j, k);
                double vzeta=dv.getdeta().getValue(i, j, k);
                
                double wxi=dw.getdxi().getValue(i, j, k);
                double weta=dw.getdeta().getValue(i, j, k);
                double wzeta=dw.getdeta().getValue(i, j, k);
                
                double txi=dT.getdxi().getValue(i, j, k);
                double teta=dT.getdeta().getValue(i, j, k);
                double tzeta=dT.getdeta().getValue(i, j, k);
                
                //Jacobians
                double xi_x=xi_xyz[0].getValue(i, j, k);
                double eta_x=xi_xyz[1].getValue(i, j, k);
                double zeta_x=xi_xyz[2].getValue(i, j, k);
                
                double xi_y=eta_xyz[0].getValue(i, j, k);
                double eta_y=eta_xyz[1].getValue(i, j, k);
                double zeta_y=eta_xyz[2].getValue(i, j, k);
                
                double xi_z=zeta_xyz[0].getValue(i, j, k);
                double eta_z=zeta_xyz[1].getValue(i, j, k);
                double zeta_z=zeta_xyz[2].getValue(i, j, k);
                
                double jac=jacobian.getValue(i, j, k);
                
                //Miu and Coeff
                double miu0=miu.getValue(i, j, k);
                double coeff0=coeff.getValue(i, j, k);
                
                //d(u,v,w,T)/d(x,y,z)
                double ux=uxi*xi_x+ueta*eta_x+uzeta*zeta_x;
                double uy=uxi*xi_y+ueta*eta_y+uzeta*zeta_y;
                double uz=uxi*xi_z+ueta*eta_z+uzeta*zeta_z;
                
                double vx=vxi*xi_x+veta*eta_x+vzeta*zeta_x;
                double vy=vxi*xi_y+veta*eta_y+vzeta*zeta_y;
                double vz=vxi*xi_z+veta*eta_z+vzeta*zeta_z;
                
                double wx=wxi*xi_x+weta*eta_x+wzeta*zeta_x;
                double wy=wxi*xi_y+weta*eta_y+wzeta*zeta_y;
                double wz=wxi*xi_z+weta*eta_z+wzeta*zeta_z;
                
                double tx=txi*xi_x+teta*eta_x+tzeta*zeta_x;
                double ty=txi*xi_y+teta*eta_y+tzeta*zeta_y;
                double tz=txi*xi_z+teta*eta_z+tzeta*zeta_z;
                
                double grad=ux+vy+wz;
                
                double tau_xx=miu0*(2.0*ux-2.0/3.0*grad);
                double tau_xy=miu0*(uy+vx-2.0/3.0*grad);
                double tau_xz=miu0*(uz+wx-2.0/3.0*grad);
                double tau_yy=miu0*(2.0*vy-2.0/3.0*grad);
                double tau_yz=miu0*(vz+wy-2.0/3.0*grad);
                double tau_zz=miu0*(2.0*wz-2.0/3.0*grad);
                
                double tau_tx=coeff0*tx;
                double tau_ty=coeff0*ty;
                double tau_tz=coeff0*tz;
                
//                tau[0].setValue(i,j, k, tau_xx);
//                tau[1].setValue(i,j, k, tau_xy);
//                tau[2].setValue(i,j, k, tau_xz);
//                tau[3].setValue(i,j, k, tau_yy);
//                tau[4].setValue(i,j, k, tau_yz);
//                tau[5].setValue(i,j, k, tau_zz);
//                tau[6].setValue(i,j, k, tau_tx);
//                tau[7].setValue(i,j, k, tau_ty);
//                tau[8].setValue(i,j, k, tau_tz);
                
                double fv[5];
                double gv[5];
                double hv[5];
                
                double fv_xi[5];
                double fv_eta[5];
                double fv_zeta[5];
                
                fv[0]=0.0;
                fv[1]=tau_xx;
                fv[2]=tau_xy;
                fv[3]=tau_xz;
                fv[4]=u*tau_xx+v*tau_xy+w*tau_xz+tau_tx;
                
                gv[0]=0.0;
                gv[1]=tau_xy;
                gv[2]=tau_yy;
                gv[3]=tau_yz;
                gv[4]=u*tau_xy+v*tau_yy+w*tau_yz+tau_ty;
                
                hv[0]=0.0;
                hv[1]=tau_xz;
                hv[2]=tau_yz;
                hv[3]=tau_zz;
                hv[4]=u*tau_xz+v*tau_yz+w*tau_zz+tau_tz;
                
                for(int iter=0;iter<5;iter++)
                {
                    fv_xi[iter]=(xi_x*fv[iter]+xi_y*gv[iter]+xi_z*hv[iter])/jac;
                    fv_eta[iter]=(eta_x*fv[iter]+eta_y*gv[iter]+eta_z*hv[iter])/jac;
                    fv_zeta[iter]=(zeta_x*fv[iter]+zeta_y*gv[iter]+zeta_z*hv[iter])/jac;
                }
                
                Flux_xi_vis[0].setValue(i, j, k,fv_xi[0]);
                Flux_xi_vis[1].setValue(i, j, k,fv_xi[1]);
                Flux_xi_vis[2].setValue(i, j, k,fv_xi[2]);
                Flux_xi_vis[3].setValue(i, j, k,fv_xi[3]);
                Flux_xi_vis[4].setValue(i, j, k,fv_xi[4]);
                
                Flux_eta_vis[0].setValue(i, j, k,fv_eta[0]);
                Flux_eta_vis[1].setValue(i, j, k,fv_eta[1]);
                Flux_eta_vis[2].setValue(i, j, k,fv_eta[2]);
                Flux_eta_vis[3].setValue(i, j, k,fv_eta[3]);
                Flux_eta_vis[4].setValue(i, j, k,fv_eta[4]);
                
                Flux_zeta_vis[0].setValue(i, j, k,fv_zeta[0]);
                Flux_zeta_vis[1].setValue(i, j, k,fv_zeta[1]);
                Flux_zeta_vis[2].setValue(i, j, k,fv_zeta[2]);
                Flux_zeta_vis[3].setValue(i, j, k,fv_zeta[3]);
                Flux_zeta_vis[4].setValue(i, j, k,fv_zeta[4]);
            }
        }
    }
    
}

void nuc3d::NaiverStokesData3d::setBoundaryViscousFlux()
{
    
}

void nuc3d::NaiverStokesData3d::setDerivativesVis()
{
    
}

void nuc3d::NaiverStokesData3d::solveRHS(PDEData3d &)
{
    
}





nuc3d::NaiverStokesData3d::~NaiverStokesData3d()
{}

