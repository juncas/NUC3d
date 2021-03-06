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
#include "gradvector.hpp"
/*************************************************************************************
 Member functions of class: NaiverStokesData3d
 **************************************************************************************/
nuc3d::NaiverStokesData3d::NaiverStokesData3d(int nx0, int ny0, int nz0, int neqs):
EulerData3D(nx0,ny0,nz0,neqs),
du(nx0,ny0,nz0),
dv(nx0,ny0,nz0),
dw(nx0,ny0,nz0),
dT(nx0,ny0,nz0),
miu(nx0,ny0,nz0),
coeff(nx0,ny0,nz0),
tau(9,Field(nx0,ny0,nz0)),
Flux_xi_vis(neqs,Field(nx0,ny0,nz0)),
Flux_eta_vis(neqs,Field(ny0,nz0,nx0)),
Flux_zeta_vis(neqs,Field(nz0,nx0,ny0)),
dfvdxi(neqs,Field(nx0,ny0,nz0)),
dgvdeta(neqs,Field(ny0,nz0,nx0)),
dhvdzeta(neqs,Field(nz0,nx0,ny0))
{
    
}

void nuc3d::NaiverStokesData3d::solve(PDEData3d &myPDE,
                                      fieldOperator3d &myOP,
                                      std::vector<bufferData> &myBf,
                                      physicsModel &myModel,
                                      MPIComunicator3d_nonblocking &myMPI,
                                      boundaryCondition &myBC)
{
    EulerData3D::solveCon2Prim(myPDE, myModel);
    EulerData3D::solveRiemann(myPDE, myModel);
    EulerData3D::setBoundaryCondition(myPDE,myModel,myBf,myBC);
    EulerData3D::solveInv(myOP,myBf,myMPI,myBC);
    NaiverStokesData3d::solveVis(myPDE,myOP,myModel,myBf,myMPI,myBC);
    NaiverStokesData3d::solveRHS(myPDE);
}

void nuc3d::NaiverStokesData3d::solveVis(PDEData3d &myPDE,
                                         fieldOperator3d &myOP,
                                         physicsModel &myModel,
                                         std::vector<bufferData> &myBf,
                                         MPIComunicator3d_nonblocking &myMPI,
                                         boundaryCondition &myBC)
{
    solveGrads(myPDE, myOP, myBf, myMPI,myBC);
    solveViscousFlux(myModel);
    setDerivativesVis(myOP,myBf,myMPI,myBC);
}

void nuc3d::NaiverStokesData3d::solveGrads(PDEData3d &myPDE,
                                           fieldOperator3d &myOP,
                                           std::vector<bufferData> &myBf,
                                           MPIComunicator3d_nonblocking &myMPI,
                                           boundaryCondition &myBC)
{
    Field &u=this->EulerData3D::W_Euler[1];
    Field &v=this->EulerData3D::W_Euler[2];
    Field &w=this->EulerData3D::W_Euler[3];
    Field &T=this->EulerData3D::W0_Euler[0];
    
    
    solve_grad(u,myOP,myBf[0],myMPI,du,0,myBC);
    
    solve_grad(v,myOP,myBf[1],myMPI,dv,1,myBC);
    
    solve_grad(w,myOP,myBf[2],myMPI,dw,2,myBC);
    
    solve_grad(T,myOP,myBf[3],myMPI,dT,3,myBC);
    
    MPI_Barrier(MPI_COMM_WORLD);
}

void nuc3d::NaiverStokesData3d::solve_grad(Field &myField,
                                           fieldOperator3d &myOP,
                                           bufferData &myBf,
                                           MPIComunicator3d_nonblocking &myMPI,
                                           gradvector &myGrad,
                                           int fdID,
                                           boundaryCondition &myBC)
{
    int typeL;
    int typeR;
    Field &dxi=myGrad.getdxi();
    Field &deta=myGrad.getdeta();
    Field &dzeta=myGrad.getdzeta();
    
    Field &fxi=myGrad.getf_xi();
    Field &feta=myGrad.getf_eta();
    Field &fzeta=myGrad.getf_zeta();
    
    typeL=myBC.getBCtype(0);
    typeR=myBC.getBCtype(1);
    
    myGrad.setGrad(myField);
    
    solveGradXi(fxi,myOP,myBf,myMPI,dxi,fdID,typeL,typeR);
    
    typeL=myBC.getBCtype(2);
    typeR=myBC.getBCtype(3);
    
    solveGradEta(feta,myOP,myBf,myMPI,deta,fdID,typeL,typeR);
    
    typeL=myBC.getBCtype(4);
    typeR=myBC.getBCtype(5);
    
    solveGradZeta(fzeta,myOP,myBf,myMPI,dzeta,fdID,typeL,typeR);
}

void nuc3d::NaiverStokesData3d::solveGradXi(Field &myField,
                                            fieldOperator3d &myOP,
                                            bufferData &myBf,
                                            MPIComunicator3d_nonblocking &myMPI,
                                            Field &dxi,
                                            int fdID,
                                            int typeL,
                                            int typeR)
{
    myMPI.bufferSendRecv(myField, myBf, 0, fdID);
    myOP.differenceInner(myField, 0, dxi);
    myMPI.waitSendRecv(myBf, 0);
    myOP.differenceBoundary(myField, myBf.BufferRecv[0], myBf.BufferRecv[1], 0, dxi,typeL,typeR);
    MPI_Barrier(MPI_COMM_WORLD);
}

void nuc3d::NaiverStokesData3d::solveGradEta(Field &myField,
                                             fieldOperator3d &myOP,
                                             bufferData &myBf,
                                             MPIComunicator3d_nonblocking &myMPI,
                                             Field &deta,
                                             int fdID,
                                             int typeL,
                                             int typeR)
{
    myMPI.bufferSendRecv(myField, myBf, 1, fdID);
    myOP.differenceInner(myField, 1, deta);
    myMPI.waitSendRecv(myBf, 1);
    myOP.differenceBoundary(myField, myBf.BufferRecv[2], myBf.BufferRecv[3], 1,deta,typeL,typeR);
    MPI_Barrier(MPI_COMM_WORLD);
}

void nuc3d::NaiverStokesData3d::solveGradZeta(Field &myField,
                                              fieldOperator3d &myOP,
                                              bufferData &myBf,
                                              MPIComunicator3d_nonblocking &myMPI,
                                              Field &dzeta,
                                              int fdID,
                                              int typeL,
                                              int typeR)
{
    myMPI.bufferSendRecv(myField, myBf, 2, fdID);
    myOP.differenceInner(myField, 2, dzeta);
    myMPI.waitSendRecv(myBf, 2);
    myOP.differenceBoundary(myField, myBf.BufferRecv[4], myBf.BufferRecv[5], 2, dzeta,typeL,typeR);
    MPI_Barrier(MPI_COMM_WORLD);
}

void nuc3d::NaiverStokesData3d::solveViscousFlux(physicsModel &myPhyMod)
{
    
    myPhyMod.getMiu(this->EulerData3D::W0_Euler[0], miu, coeff);
    
    double *pu=this->EulerData3D::W_Euler[1].getDataPtr();
    double *pv=this->EulerData3D::W_Euler[2].getDataPtr();
    double *pw=this->EulerData3D::W_Euler[3].getDataPtr();
    
    double *uxi=du.getdxi().getDataPtr();
    double *vxi=dv.getdxi().getDataPtr();
    double *wxi=dw.getdxi().getDataPtr();
    double *txi=dT.getdxi().getDataPtr();
    
    double *ueta=du.getdeta().getDataPtr();
    double *veta=dv.getdeta().getDataPtr();
    double *weta=dw.getdeta().getDataPtr();
    double *teta=dT.getdeta().getDataPtr();
    
    double *uzeta=du.getdzeta().getDataPtr();
    double *vzeta=dv.getdzeta().getDataPtr();
    double *wzeta=dw.getdzeta().getDataPtr();
    double *tzeta=dT.getdzeta().getDataPtr();
    
    double *pjac=jacobian.getDataPtr();
    double *pxi_x=xi_xyz[0].getDataPtr();
    double *pxi_y=xi_xyz[1].getDataPtr();
    double *pxi_z=xi_xyz[2].getDataPtr();
    
    double *peta_x=eta_xyz[0].getDataPtr();
    double *peta_y=eta_xyz[1].getDataPtr();
    double *peta_z=eta_xyz[2].getDataPtr();
    
    double *pzeta_x=zeta_xyz[0].getDataPtr();
    double *pzeta_y=zeta_xyz[1].getDataPtr();
    double *pzeta_z=zeta_xyz[2].getDataPtr();
    
    double *pMiu=miu.getDataPtr();
    double *pCoeff=coeff.getDataPtr();
    
    double *flux_xi[5];
    double *flux_eta[5];
    double *flux_zeta[5];
    
    double fv[5];
    double gv[5];
    double hv[5];
    
    flux_xi[0]=Flux_xi_vis[0].getDataPtr();
    flux_xi[1]=Flux_xi_vis[1].getDataPtr();
    flux_xi[2]=Flux_xi_vis[2].getDataPtr();
    flux_xi[3]=Flux_xi_vis[3].getDataPtr();
    flux_xi[4]=Flux_xi_vis[4].getDataPtr();
    
    flux_eta[0]=Flux_eta_vis[0].getDataPtr();
    flux_eta[1]=Flux_eta_vis[1].getDataPtr();
    flux_eta[2]=Flux_eta_vis[2].getDataPtr();
    flux_eta[3]=Flux_eta_vis[3].getDataPtr();
    flux_eta[4]=Flux_eta_vis[4].getDataPtr();
    
    flux_zeta[0]=Flux_zeta_vis[0].getDataPtr();
    flux_zeta[1]=Flux_zeta_vis[1].getDataPtr();
    flux_zeta[2]=Flux_zeta_vis[2].getDataPtr();
    flux_zeta[3]=Flux_zeta_vis[3].getDataPtr();
    flux_zeta[4]=Flux_zeta_vis[4].getDataPtr();
    
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                int idx_eta=ny*nz*i+ny*k+j;
                int idx_zeta=nz*nx*j+nz*i+k;
                
                double jac=pjac[idx_xi];
                double xi_x=pxi_x[idx_xi];
                double xi_y=pxi_y[idx_xi];
                double xi_z=pxi_z[idx_xi];
                
                double eta_x=peta_x[idx_xi];
                double eta_y=peta_y[idx_xi];
                double eta_z=peta_z[idx_xi];
                
                double zeta_x=pzeta_x[idx_xi];
                double zeta_y=pzeta_y[idx_xi];
                double zeta_z=pzeta_z[idx_xi];
                
                //d(u,v,w,T)/d(x,y,z)
                double ux=uxi[idx_xi]*xi_x+ueta[idx_eta]*eta_x+uzeta[idx_zeta]*zeta_x;
                double uy=uxi[idx_xi]*xi_y+ueta[idx_eta]*eta_y+uzeta[idx_zeta]*zeta_y;
                double uz=uxi[idx_xi]*xi_z+ueta[idx_eta]*eta_z+uzeta[idx_zeta]*zeta_z;
                
                double vx=vxi[idx_xi]*xi_x+veta[idx_eta]*eta_x+vzeta[idx_zeta]*zeta_x;
                double vy=vxi[idx_xi]*xi_y+veta[idx_eta]*eta_y+vzeta[idx_zeta]*zeta_y;
                double vz=vxi[idx_xi]*xi_z+veta[idx_eta]*eta_z+vzeta[idx_zeta]*zeta_z;
                
                double wx=wxi[idx_xi]*xi_x+weta[idx_eta]*eta_x+wzeta[idx_zeta]*zeta_x;
                double wy=wxi[idx_xi]*xi_y+weta[idx_eta]*eta_y+wzeta[idx_zeta]*zeta_y;
                double wz=wxi[idx_xi]*xi_z+weta[idx_eta]*eta_z+wzeta[idx_zeta]*zeta_z;
                
                double tx=txi[idx_xi]*xi_x+teta[idx_eta]*eta_x+tzeta[idx_zeta]*zeta_x;
                double ty=txi[idx_xi]*xi_y+teta[idx_eta]*eta_y+tzeta[idx_zeta]*zeta_y;
                double tz=txi[idx_xi]*xi_z+teta[idx_eta]*eta_z+tzeta[idx_zeta]*zeta_z;
                
                double grad=ux+vy+wz;
                double miu0=pMiu[idx_xi];
                
                double tau_xx=miu0*(2.0*ux-2.0/3.0*grad);
                double tau_xy=miu0*(uy+vx-2.0/3.0*grad);
                double tau_xz=miu0*(uz+wx-2.0/3.0*grad);
                double tau_yy=miu0*(2.0*vy-2.0/3.0*grad);
                double tau_yz=miu0*(vz+wy-2.0/3.0*grad);
                double tau_zz=miu0*(2.0*wz-2.0/3.0*grad);
                
                double coeff0=pCoeff[idx_xi];
                double tau_tx=coeff0*tx;
                double tau_ty=coeff0*ty;
                double tau_tz=coeff0*tz;                
                
                double u=pu[idx_xi];
                double v=pv[idx_xi];
                double w=pw[idx_xi];
                
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
                    flux_xi[iter][idx_xi]=(xi_x*fv[iter]+xi_y*gv[iter]+xi_z*hv[iter])/jac;
                    flux_eta[iter][idx_eta]=(eta_x*fv[iter]+eta_y*gv[iter]+eta_z*hv[iter])/jac;
                    flux_zeta[iter][idx_zeta]=(zeta_x*fv[iter]+zeta_y*gv[iter]+zeta_z*hv[iter])/jac;
                }
                
            }
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

}

void nuc3d::NaiverStokesData3d::setDerivativesVis(fieldOperator3d &myOP,
                                                  std::vector<bufferData> &myBf,
                                                  MPIComunicator3d_nonblocking &myMPI,
                                                  boundaryCondition &myBC)
{
    setDerivativeXi(myOP,myBf,myMPI,myBC);
    setDerivativeEta(myOP,myBf,myMPI,myBC);
    setDerivativeZeta(myOP,myBf,myMPI,myBC);
    MPI_Barrier(MPI_COMM_WORLD);
}

void nuc3d::NaiverStokesData3d::setDerivativeXi(fieldOperator3d &myOP,
                                                std::vector<bufferData> &myBf,
                                                MPIComunicator3d_nonblocking &myMPI,
                                                boundaryCondition &myBC)
{
    auto beg=Flux_xi_vis.begin();
    auto end=Flux_xi_vis.end();
    int typeL=myBC.getBCtype(0);
    int typeR=myBC.getBCtype(1);
    for(auto iter=beg;iter!=end;iter++)
    {
        myMPI.bufferSendRecv(*iter, myBf[iter-beg], 0, static_cast<int>(iter-beg));
        myOP.differenceInner(*iter, 0, dfvdxi[iter-beg]);
        myMPI.waitSendRecv(myBf[iter-beg], 0);
        myOP.differenceBoundary(*iter, myBf[iter-beg].BufferRecv[0], myBf[iter-beg].BufferRecv[1], 0, dfvdxi[iter-beg],typeL,typeR);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
}

void nuc3d::NaiverStokesData3d::setDerivativeEta(fieldOperator3d &myOP,
                                                 std::vector<bufferData> &myBf,
                                                 MPIComunicator3d_nonblocking &myMPI,
                                                 boundaryCondition &myBC)
{
    auto beg=Flux_eta_vis.begin();
    auto end=Flux_eta_vis.end();
    int typeL=myBC.getBCtype(2);
    int typeR=myBC.getBCtype(3);
    for(auto iter=beg;iter!=end;iter++)
    {
        myMPI.bufferSendRecv(*iter, myBf[iter-beg], 1, static_cast<int>(iter-beg));
        myOP.differenceInner(*iter, 1, dgvdeta[iter-beg]);
        myMPI.waitSendRecv(myBf[iter-beg], 1);
        myOP.differenceBoundary(*iter, myBf[iter-beg].BufferRecv[2], myBf[iter-beg].BufferRecv[3], 1, dgvdeta[iter-beg],typeL,typeR);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
}

void nuc3d::NaiverStokesData3d::setDerivativeZeta(fieldOperator3d &myOP,
                                                  std::vector<bufferData> &myBf,
                                                  MPIComunicator3d_nonblocking &myMPI,
                                                  boundaryCondition &myBC)
{
    auto beg=Flux_zeta_vis.begin();
    auto end=Flux_zeta_vis.end();
    int typeL=myBC.getBCtype(4);
    int typeR=myBC.getBCtype(5);
    for(auto iter=beg;iter!=end;iter++)
    {
        myMPI.bufferSendRecv(*iter, myBf[iter-beg], 2, static_cast<int>(iter-beg));
        myOP.differenceInner(*iter, 2, dhvdzeta[iter-beg]);
        myMPI.waitSendRecv(myBf[iter-beg], 2);
        myOP.differenceBoundary(*iter, myBf[iter-beg].BufferRecv[4], myBf[iter-beg].BufferRecv[5], 2, dhvdzeta[iter-beg],typeL,typeR);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);    
}


void nuc3d::NaiverStokesData3d::solveRHS(PDEData3d &myPDE)
{
    VectorField &df_v=this->EulerData3D::dfdxi;
    VectorField &dg_v=this->EulerData3D::dgdeta;
    VectorField &dh_v=this->EulerData3D::dhdzeta;
    
    VectorField &dfv_v=this->dfvdxi;
    VectorField &dgv_v=this->dgvdeta;
    VectorField &dhv_v=this->dhvdzeta;
    
    auto beg=myPDE.getRHS().begin();
    auto end=myPDE.getRHS().end();
    
    for (auto iter=beg ; iter!=end ; iter++)
    {
        double *df=df_v[iter-beg].getDataPtr();
        double *dg=dg_v[iter-beg].getDataPtr();
        double *dh=dh_v[iter-beg].getDataPtr();
        
        double *dfv=dfv_v[iter-beg].getDataPtr();
        double *dgv=dgv_v[iter-beg].getDataPtr();
        double *dhv=dhv_v[iter-beg].getDataPtr();
        
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
        
        for (int k=0; k<nz; k++)
        {
            for (int j=0; j<ny; j++)
            {
                for (int i=0; i<nx; i++)
                {
                    int idx_xi=nx*ny*k+nx*j+i;
                    
                    rhs[idx_xi]-=dfv[idx_xi];
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
                    
                    rhs[idx_xi]-=dgv[idx_eta];
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
                    
                    rhs[idx_xi]-=dhv[idx_zeta];
                }
            }
        }
        
    }
    
    this->EulerData3D::getDt();
    myPDE.setDt(dt);
    
}

nuc3d::NaiverStokesData3d::~NaiverStokesData3d()
{}

