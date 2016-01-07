//
//  postproc.cpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/4.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//

#include "postproc.hpp"
#include "physicsModel.h"
#include "MPICommunicator.h"
#include "IOcontroller.h"
#include "boundaryConditions.hpp"
#include "fieldOperator.h"

nuc3d::postproc::postproc(int nx0,int ny0,int nz0):
nx(nx0),
ny(ny0),
nz(nz0),
Q(nx0,ny0,nz0),
omega0(nx0,ny0,nz0),
omega1(nx0,ny0,nz0),
omega2(nx0,ny0,nz0),
q_xi(nx0,ny0,nz0),
q_eta(ny0,nz0,nx0),
q_zeta(nz0,nx0,ny0),
u_xi(nx0,ny0,nz0),
u_eta(ny0,nz0,nx0),
u_zeta(nz0,nx0,ny0),
v_xi(nx0,ny0,nz0),
v_eta(ny0,nz0,nx0),
v_zeta(nz0,nx0,ny0),
w_xi(nx0,ny0,nz0),
w_eta(ny0,nz0,nx0),
w_zeta(nz0,nx0,ny0),
u_x(nx0,ny0,nz0),
u_y(nx0,ny0,nz0),
u_z(nx0,ny0,nz0),
v_x(nx0,ny0,nz0),
v_y(nx0,ny0,nz0),
v_z(nx0,ny0,nz0),
w_x(nx0,ny0,nz0),
w_y(nx0,ny0,nz0),
w_z(nx0,ny0,nz0),
enstrophy(0.0),
postSize(4)
{
    
}

nuc3d::postproc::~postproc()
{
    
}

void nuc3d::postproc::solvePost(VectorField &prims,
                                VectorField &acous,
                                VectorField &xyz,
                                VectorField &xi_xyz,
                                VectorField &eta_xyz,
                                VectorField &zeta_xyz,
                                fieldOperator3d &myOP,
                                VectorBuffer &myBf,
                                MPIComunicator3d_nonblocking &myMPI,
                                boundaryCondition &myBC,
                                IOController &myIO,
                                int istep,
                                double time)
{
    solveQ(prims, xi_xyz, eta_xyz, zeta_xyz, myOP, myBf, myMPI, myBC);
    
    MPI_Allreduce(&enstrophy, &enstrophy_glb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    MPI_Allreduce(&kinetic, &kinetic_glb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    if(0==myMPI.getMyId())
    {
        std::ofstream myIOfile;
        myIOfile.open("enstrophy.dat",std::ios::out|std::ios::app);
        myIOfile<<std::setprecision(12)
        <<time
        <<" "
        <<enstrophy_glb
        <<" "
        <<kinetic_glb<<"\n";
        myIOfile.close();
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
}

void nuc3d::postproc::OutputPost(VectorField &prims,
                                 VectorField &acous,
                                 VectorField &xyz,
                                 VectorField &xi_xyz,
                                 VectorField &eta_xyz,
                                 VectorField &zeta_xyz,
                                 fieldOperator3d &myOP,
                                 VectorBuffer &myBf,
                                 MPIComunicator3d_nonblocking &myMPI,
                                 boundaryCondition &myBC,
                                 IOController &myIO,
                                 int istep,
                                 double time)
{
    std::string forename_flow = ("flowData/Q_step_");
    std::string step;
    std::string mid("_id_");
    std::string id;
    std::string tailname = (".dat");
    int myID=myMPI.getMyId();
    
    std::stringstream ss_step,ss_id;
    ss_step << istep;
    ss_step >> step;
    
    ss_id<<myID;
    ss_id>>id;
    
    std::string filename_flow = forename_flow + step + mid + id + tailname;
    
    std::ofstream myIOfile;
    
    myIOfile.open(filename_flow);
    
    std::string TECplotHeader[2]={"title=time_",
        "variables=x,y,z"};
    
    myIOfile<<TECplotHeader[0]<<time<<"\n"
    <<TECplotHeader[1];
    
    for(int i=0;i<postSize;i++)
    {
        std::string head("Val_");
        std::string temp;
        std::stringstream sstemp;
        
        sstemp<<i;
        sstemp>>temp;
        myIOfile<<","<<head+temp;
    }
    
    myIOfile<<"\n Zone I = "<<nx+1<<", J= "<<ny+1<<", K="<<nz+1
    <<"\n DATAPACKING=BLOCK, VARLOCATION=(["<<xyz.size()+1<<"-"
    <<xyz.size()+postSize
    <<"]=CELLCENTERED)\n";
    
    for(auto iter=xyz.begin();iter!=xyz.end();iter++)
    {
        writeField(myIOfile, *iter);
    }
    
    writeField(myIOfile, Q);
    writeField(myIOfile, omega0);
    writeField(myIOfile, omega1);
    writeField(myIOfile, omega2);
    
    /* ADD arbitrary number of post values
     writeField(myIOfile, omega2);
     */
    
    
    myIOfile.close();
}
void nuc3d::postproc::solveQ(VectorField &prims,
                             VectorField &xi_xyz,
                             VectorField &eta_xyz,
                             VectorField &zeta_xyz,
                             fieldOperator3d &myOP,
                             VectorBuffer &myBf,
                             MPIComunicator3d_nonblocking &myMPI,
                             boundaryCondition &myBC)
{
    Field &rho=prims[0];
    Field &u=prims[1];
    Field &v=prims[2];
    Field &w=prims[3];
    
    solveGrad(u, myOP, myBf[0], myMPI,myBC,0,u_xi,u_eta,u_zeta);
    solveGrad(v, myOP, myBf[1], myMPI,myBC,1,v_xi,v_eta,v_zeta);
    solveGrad(w, myOP, myBf[2], myMPI,myBC,2,w_xi,w_eta,w_zeta);
    solveGrad_xyz(xi_xyz, eta_xyz, zeta_xyz);
    
    
    int nx=Q.getSizeX();
    int ny=Q.getSizeY();
    int nz=Q.getSizeZ();
    
    double *prho=rho.getDataPtr();
    double *pu=u.getDataPtr();
    double *pv=v.getDataPtr();
    double *pw=w.getDataPtr();
    
    double *ux=u_x.getDataPtr();
    double *uy=u_y.getDataPtr();
    double *uz=u_z.getDataPtr();
    
    double *vx=v_x.getDataPtr();
    double *vy=v_y.getDataPtr();
    double *vz=v_z.getDataPtr();
    
    double *wx=w_x.getDataPtr();
    double *wy=w_y.getDataPtr();
    double *wz=w_z.getDataPtr();
    
    double *pOmega0=omega0.getDataPtr();
    double *pOmega1=omega1.getDataPtr();
    double *pOmega2=omega2.getDataPtr();
    double *pQ=Q.getDataPtr();
    
    enstrophy=0.0;
    kinetic=0.0;
    
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                double w0,w1,w2;
                double q0;
                
                w0=wy[idx_xi]-vz[idx_xi];
                w1=uz[idx_xi]-wx[idx_xi];
                w2=vx[idx_xi]-uy[idx_xi];
                
                q0=ux[idx_xi]*vy[idx_xi]
                +ux[idx_xi]*wz[idx_xi]
                +vy[idx_xi]*wz[idx_xi]
                -uy[idx_xi]*vx[idx_xi]
                -uz[idx_xi]*wx[idx_xi]
                -vz[idx_xi]*wy[idx_xi];
                
                pOmega0[idx_xi]=w0;
                pOmega1[idx_xi]=w1;
                pOmega2[idx_xi]=w2;
                pQ[idx_xi]=q0;
                enstrophy+=0.5*(w0*w0+w1*w1+w2*w2);
                kinetic+=0.5*prho[idx_xi]*(pu[idx_xi]*pu[idx_xi]
                                           +pv[idx_xi]*pv[idx_xi]
                                           +pw[idx_xi]*pw[idx_xi]);
            }
        }
    }
    
}

void nuc3d::postproc::solveGrad(Field &myField,
                                fieldOperator3d &myOP,
                                bufferData &myBf,
                                MPIComunicator3d_nonblocking &myMPI,
                                boundaryCondition &myBC,
                                int fdID,
                                Field &dfdxi,
                                Field &dfdeta,
                                Field &dfdzeta)
{
    int typeL;
    int typeR;
    
    setField(myField);
    
    typeL=myBC.getBCtype(0);
    typeR=myBC.getBCtype(1);
    solveGrad_df(q_xi, myOP, myBf, myMPI, dfdxi, 0, fdID, typeL, typeR);
    
    typeL=myBC.getBCtype(2);
    typeR=myBC.getBCtype(3);
    solveGrad_df(q_eta, myOP, myBf, myMPI, dfdeta, 1, fdID, typeL, typeR);
    
    typeL=myBC.getBCtype(4);
    typeR=myBC.getBCtype(5);
    solveGrad_df(q_zeta, myOP, myBf, myMPI, dfdzeta, 2, fdID, typeL, typeR);
    
    
}

void nuc3d::postproc::solveGrad_xyz(VectorField &xi_xyz,
                                    VectorField &eta_xyz,
                                    VectorField &zeta_xyz)
{
    int nx=Q.getSizeX();
    int ny=Q.getSizeY();
    int nz=Q.getSizeZ();
    
    double *pxi_x=xi_xyz[0].getDataPtr();
    double *pxi_y=xi_xyz[1].getDataPtr();
    double *pxi_z=xi_xyz[2].getDataPtr();
    
    double *peta_x=eta_xyz[0].getDataPtr();
    double *peta_y=eta_xyz[1].getDataPtr();
    double *peta_z=eta_xyz[2].getDataPtr();
    
    double *pzeta_x=zeta_xyz[0].getDataPtr();
    double *pzeta_y=zeta_xyz[1].getDataPtr();
    double *pzeta_z=zeta_xyz[2].getDataPtr();
    
    double *uxi=u_xi.getDataPtr();
    double *ueta=u_eta.getDataPtr();
    double *uzeta=u_zeta.getDataPtr();
    
    double *vxi=v_xi.getDataPtr();
    double *veta=v_eta.getDataPtr();
    double *vzeta=v_zeta.getDataPtr();
    
    double *wxi=w_xi.getDataPtr();
    double *weta=w_eta.getDataPtr();
    double *wzeta=w_zeta.getDataPtr();
    
    double *ux=u_x.getDataPtr();
    double *uy=u_y.getDataPtr();
    double *uz=u_z.getDataPtr();
    
    double *vx=v_x.getDataPtr();
    double *vy=v_y.getDataPtr();
    double *vz=v_z.getDataPtr();
    
    double *wx=w_x.getDataPtr();
    double *wy=w_y.getDataPtr();
    double *wz=w_z.getDataPtr();
    
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                int idx_eta=ny*nz*i+ny*k+j;
                int idx_zeta=nz*nx*j+nz*i+k;
                
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
                ux[idx_xi]=uxi[idx_xi]*xi_x+ueta[idx_eta]*eta_x+uzeta[idx_zeta]*zeta_x;
                uy[idx_xi]=uxi[idx_xi]*xi_y+ueta[idx_eta]*eta_y+uzeta[idx_zeta]*zeta_y;
                uz[idx_xi]=uxi[idx_xi]*xi_z+ueta[idx_eta]*eta_z+uzeta[idx_zeta]*zeta_z;
                
                vx[idx_xi]=vxi[idx_xi]*xi_x+veta[idx_eta]*eta_x+vzeta[idx_zeta]*zeta_x;
                vy[idx_xi]=vxi[idx_xi]*xi_y+veta[idx_eta]*eta_y+vzeta[idx_zeta]*zeta_y;
                vz[idx_xi]=vxi[idx_xi]*xi_z+veta[idx_eta]*eta_z+vzeta[idx_zeta]*zeta_z;
                
                wx[idx_xi]=wxi[idx_xi]*xi_x+weta[idx_eta]*eta_x+wzeta[idx_zeta]*zeta_x;
                wy[idx_xi]=wxi[idx_xi]*xi_y+weta[idx_eta]*eta_y+wzeta[idx_zeta]*zeta_y;
                wz[idx_xi]=wxi[idx_xi]*xi_z+weta[idx_eta]*eta_z+wzeta[idx_zeta]*zeta_z;
            }
        }
    }
    
}

void nuc3d::postproc::setField(const Field &myField)
{
    
    double *pField=myField.getDataPtr();
    double *pf_xi=q_xi.getDataPtr();
    double *pf_eta=q_eta.getDataPtr();
    double *pf_zeta=q_zeta.getDataPtr();
    
    int nx0=myField.getSizeX();
    int ny0=myField.getSizeY();
    int nz0=myField.getSizeZ();
    
    for (int k=0; k<nz0; k++)
    {
        for (int j=0; j<ny0; j++)
        {
            for (int i=0; i<nx0; i++)
            {
                int idx=nx0*ny0*k+nx0*j+i;
                pf_xi[idx]=pField[idx];
            }
        }
    }
    
    for (int k=0; k<nz0; k++)
    {
        for (int j=0; j<ny0; j++)
        {
            for (int i=0; i<nx0; i++)
            {
                int idx=nx0*ny0*k+nx0*j+i;
                int idx_eta=ny0*nz0*i+ny0*k+j;
                
                pf_eta[idx_eta]=pField[idx];
            }
        }
    }
    
    for (int k=0; k<nz0; k++)
    {
        for (int j=0; j<ny0; j++)
        {
            for (int i=0; i<nx0; i++)
            {
                int idx=nx0*ny0*k+nx0*j+i;
                int idx_zeta=nz0*nx0*j+nz0*i+k;
                
                pf_zeta[idx_zeta]=pField[idx];
            }
        }
    }
    
}

void nuc3d::postproc::solveGrad_df(Field &myField,
                                   fieldOperator3d &myOP,
                                   bufferData &myBf,
                                   MPIComunicator3d_nonblocking &myMPI,
                                   Field &df,
                                   int dir,
                                   int fdID,
                                   int typeL,
                                   int typeR)
{
    myMPI.bufferSendRecv(myField, myBf, dir, fdID);
    myOP.differenceInner(myField, dir, df);
    myMPI.waitSendRecv(myBf, dir);
    myOP.differenceBoundary(myField, myBf.BufferRecv[2*dir], myBf.BufferRecv[2*dir+1], dir, df,typeL,typeR);
    MPI_Barrier(MPI_COMM_WORLD);
}


void nuc3d::postproc::writeField(std::ofstream &myFile, Field &myField)
{
    int nx0=myField.getSizeX();
    int ny0=myField.getSizeY();
    int nz0=myField.getSizeZ();
    
    for(int k=0;k<nz0;k++)
    {
        for(int j=0;j<ny0;j++)
        {
            for(int i=0;i<nx0;i++)
            {
                double value=myField.getValue(i,j,k);
                myFile<<std::setprecision(12)<<value<<"\n";
            }
        }
    }
}

