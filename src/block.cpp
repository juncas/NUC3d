//
//  block.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/10/20.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#include "block.h"
#include "euler3d.h"
#include "eulerReactive3d.h"
#include "NaiverStokes3d.h"
#include "NaiverStokesReactive3d.h"
#include "physicsModel.h"
#include "MPICommunicator.h"
#include "IOcontroller.h"
#include "boundaryConditions.hpp"
#include "postproc.hpp"

#define TILE_SIZE 4;

nuc3d::block::block():
time(0.0),
dt(0.025),
istep(0),
RES(1.0)
{}

nuc3d::block::~block()
{}

void nuc3d::block::initial(fieldOperator3d &myOP,
                           physicsModel &myPhyMod,
                           MPIComunicator3d_nonblocking &myMPI,
                           boundaryCondition &myBC,
                           IOController &myIO)

{
    std::stringstream ss;
    int ProcId = myMPI.getMyId();
    int nx0, ny0, nz0;
    
    std::string forename_mesh = ("mesh/mesh_");
    std::string forename_flow = ("flow_");
    std::string midname;
    std::string tailname = (".x");
    
    ss << ProcId;
    ss >> midname;
    
    std::string filename_mesh = forename_mesh + midname + tailname;
    std::string filename_flow = forename_flow + midname + tailname;
    
    std::ifstream myFile;
    myFile.open(filename_mesh);
    if (myFile)
    {
        myFile >> nx0 >> ny0 >> nz0>> bfsize;
        nx=nx0-1;
        ny=ny0-1;
        nz=nz0-1;
        for(int i=0;i<3;i++)
        {
            xyz.push_back(Field(nx0,ny0,nz0));
            xyz_center.push_back(Field(nx,ny,nz));
        }
        if(myMPI.getMyId()==0) std::cout<<"Start reading mesh data..."<<std::endl;
        for(auto iter=xyz.begin();iter!=xyz.end();iter++)
            readField(myFile,*iter);
        if(myMPI.getMyId()==0) std::cout<<"Mesh data has been read!"<<std::endl;
        
        if(myMPI.getMyId()==0) std::cout<<"Start calculating mesh data..."<<std::endl;
        getXYZ_center();
        if(myMPI.getMyId()==0) std::cout<<"Center data has been calculated..."<<std::endl;
        //writeField(myFile_o,*iter);
    }
    else
    {
        std::cout<<"File \'"<<filename_mesh<<"\' does not exist!"<<std::endl;
        exit(-1);
    }
    myFile.close();
    
    if(myMPI.getMyId()==0) std::cout<<"Jacobians has been calculated!"<<std::endl;
    initialData(nx, ny, nz, myPhyMod);
    myBC.initialBC(mybuffer,myMPI);
    
    getJacobians(myOP,myMPI,myBC);
    
    outputGEO_tecplot(myMPI.getMyId());
    //    for(auto iter=myFluxes->getZeta_xyz().begin();iter!=myFluxes->getZeta_xyz().end();iter++)
    //      writeField(myFile_o, *iter);
    
    
    if(0==(myIO.getStep("startStep")))
        initialQ(myIO,myPhyMod);
    else
        inputQ_binary(myMPI.getMyId(),myIO.getStep("startStep"));
    
    if(myMPI.getMyId()==0) std::cout<<"Flow field has been initialized!"<<std::endl;
    
    
}

void nuc3d::block::writeField(std::ofstream &myFile, nuc3d::Field &myField)
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


void nuc3d::block::readField(std::ifstream &myFile, nuc3d::Field &myField)
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
                double value;
                if(!(myFile>>value))
                {
                    std::cout<<"Error:End of File at "
                    <<"i= "<<i<<", j="<<j<<", k= "<<k
                    <<std::endl;
                    exit(-1);
                }
                myField.setValue(i,j,k,value);
            }
        }
    }
}

void nuc3d::block::getXYZ_center()
{
    auto beg=xyz.begin();
    auto end=xyz.end();
    
    for(auto iter=beg;iter!=end;iter++)
        interpolation_lag(*iter, xyz_center[iter-beg]);
    
}

void nuc3d::block::interpolation_lag(const nuc3d::Field &input,nuc3d::Field &output)
{
    int Dim_node_xi=input.getSizeX();
    int Dim_node_eta=input.getSizeY();
    int Dim_node_zeta=input.getSizeZ();
    
    int tileDim_xi=(Dim_node_xi-1)/TILE_SIZE;
    int tileDim_eta=(Dim_node_eta-1)/TILE_SIZE;
    int tileDim_zeta=(Dim_node_zeta-1)/TILE_SIZE;
    
    for(int k_tile=0;k_tile!=tileDim_zeta;k_tile++)
    {
        for(int j_tile=0;j_tile!=tileDim_eta;j_tile++)
        {
            for(int i_tile=0;i_tile!=tileDim_xi;i_tile++)
            {
                double value=0.0;
                int ibeg_center=i_tile*TILE_SIZE;
                int jbeg_center=j_tile*TILE_SIZE;
                int kbeg_center=k_tile*TILE_SIZE;
                
                int iend_center=(i_tile+1)*TILE_SIZE;
                int jend_center=(j_tile+1)*TILE_SIZE;
                int kend_center=(k_tile+1)*TILE_SIZE;
                
                int ibeg_node=i_tile*TILE_SIZE;
                int jbeg_node=j_tile*TILE_SIZE;
                int kbeg_node=k_tile*TILE_SIZE;
                
                int iend_node=(i_tile+1)*TILE_SIZE;
                int jend_node=(j_tile+1)*TILE_SIZE;
                int kend_node=(k_tile+1)*TILE_SIZE;
                
                for(int k=kbeg_center;k<kend_center;k++)
                {
                    for(int j=jbeg_center;j<jend_center;j++)
                    {
                        for(int i=ibeg_center;i<iend_center;i++)
                        {
                            value=interpolation_lag_center(ibeg_node, iend_node,
                                                           jbeg_node, jend_node,
                                                           kbeg_node, kend_node,
                                                           input,
                                                           i-ibeg_center,j-jbeg_center, k-kbeg_center);
                            output.setValue(i, j, k, value);
                            
                        }
                    }
                }
                
            }
        }
    }
}

double nuc3d::block::interpolation_lag_center(const int ibeg,const int iend,
                                              const int jbeg,const int jend,
                                              const int kbeg,const int kend,
                                              const Field &myfield,
                                              int idx_xi,int idx_eta,int idx_zeta)
{
    double value_zeta=0.0;
    for(int k=kbeg;k<=kend;k++)
    {
        double value_eta=0.0;
        for(int j=jbeg;j<=jend;j++)
        {
            double value_xi=0.0;
            for(int i=ibeg;i<=iend;i++)
            {
                value_xi+=lag_coeff[idx_xi][i-ibeg]*myfield.getValue(i, j, k);
            }
            value_eta+=lag_coeff[idx_eta][j-jbeg]*value_xi;
        }
        value_zeta+=lag_coeff[idx_zeta][k-kbeg]*value_eta;
    }
    
    return value_zeta;
}

void nuc3d::block::getJacobians(fieldOperator3d &myOP,
                                MPIComunicator3d_nonblocking &myMPI,
                                boundaryCondition &myBC)
{
    VectorField &xi_xyz=myFluxes->getXi_xyz(); // xyz_Xi
    VectorField &eta_xyz=myFluxes->getEta_xyz();// xyz_Eta
    VectorField &zeta_xyz=myFluxes->getZeta_xyz();// xyz_Zeta
    
    VectorField &dx=myFluxes->getDx(); // xyz_Xi
    VectorField &dy=myFluxes->getDy();// xyz_Eta
    VectorField &dz=myFluxes->getDz();// xyz_Zeta
    
    Field &jac=myFluxes->getJac();
    
    //std::cout<<"Start calculate xyz_xi, xyz_eta, xyz_zta ..."<<std::endl;
    
    gradvector gradx(nx,ny,nz);
    gradvector grady(nx,ny,nz);
    gradvector gradz(nx,ny,nz);
    
    solve_grad(xyz[0], myOP, mybuffer[0], myMPI, gradx, 0, myBC);
    solve_grad(xyz[1], myOP, mybuffer[1], myMPI, grady, 1, myBC);
    solve_grad(xyz[2], myOP, mybuffer[2], myMPI, gradz, 2, myBC);
    
    //std::cout<<std::endl;
    
    //std::cout<<"xx_xyz has been calculated ..."<<std::endl;
    //std::cout<<"Calculating Jacobians calculated ..."<<std::endl;
    
    double *jacob0=jac.getDataPtr();
    
    double *xi_x0=xi_xyz[0].getDataPtr();
    double *xi_y0=xi_xyz[1].getDataPtr();
    double *xi_z0=xi_xyz[2].getDataPtr();
    
    double *eta_x0=eta_xyz[0].getDataPtr();
    double *eta_y0=eta_xyz[1].getDataPtr();
    double *eta_z0=eta_xyz[2].getDataPtr();
    
    double *zeta_x0=zeta_xyz[0].getDataPtr();
    double *zeta_y0=zeta_xyz[1].getDataPtr();
    double *zeta_z0=zeta_xyz[2].getDataPtr();
    
    double *dx0=dx[0].getDataPtr();
    double *dx1=dx[1].getDataPtr();
    double *dx2=dx[2].getDataPtr();
    
    double *dy0=dy[0].getDataPtr();
    double *dy1=dy[1].getDataPtr();
    double *dy2=dy[2].getDataPtr();
    
    double *dz0=dz[0].getDataPtr();
    double *dz1=dz[1].getDataPtr();
    double *dz2=dz[2].getDataPtr();
    

    double *x_xi0=gradx.getdxi().getDataPtr();
    double *y_xi0=grady.getdxi().getDataPtr();
    double *z_xi0=gradz.getdxi().getDataPtr();
    
    double *x_eta0=gradx.getdeta().getDataPtr();
    double *y_eta0=grady.getdeta().getDataPtr();
    double *z_eta0=gradz.getdeta().getDataPtr();
    
    double *x_zeta0=gradx.getdzeta().getDataPtr();
    double *y_zeta0=grady.getdzeta().getDataPtr();
    double *z_zeta0=gradz.getdzeta().getDataPtr();
    
    for(int k=0;k<nx;k++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int i=0;i<nx;i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                int idx_eta=ny*nz*i+ny*k+j;
                int idx_zeta=nz*nx*j+nz*i+k;
                
                double x_xi=x_xi0[idx_xi];
                double y_xi=y_xi0[idx_xi];
                double z_xi=z_xi0[idx_xi];
                
                double x_eta=x_eta0[idx_eta];
                double y_eta=y_eta0[idx_eta];
                double z_eta=z_eta0[idx_eta];
                
                double x_zeta=x_zeta0[idx_zeta];
                double y_zeta=y_zeta0[idx_zeta];
                double z_zeta=z_zeta0[idx_zeta];
                
                dx0[idx_xi]=x_xi;
                dx1[idx_xi]=x_eta;
                dx2[idx_xi]=x_zeta;
                
                dy0[idx_xi]=y_xi;
                dy1[idx_xi]=y_eta;
                dy2[idx_xi]=y_zeta;
                
                dz0[idx_xi]=z_xi;
                dz1[idx_xi]=z_eta;
                dz2[idx_xi]=z_zeta;
                
                double jacob=1.0/(x_xi*(y_eta*z_zeta-y_zeta*z_eta)
                                  +x_eta*(y_zeta*z_xi-y_xi*z_zeta)
                                  +x_zeta*(y_xi*z_eta-y_eta*z_xi));
                
                if(jacob<0.0)
                {
                    std::cout<<"Jacob < 0.0 at "<<i<<", "<<j<<", "<<k<<std::endl;
                    exit(-1);
                }

                jacob0[idx_xi]=jacob;
                xi_x0[idx_xi]=(y_eta*z_zeta-y_zeta*z_eta)*jacob0[idx_xi];
                xi_y0[idx_xi]=(x_zeta*z_eta-x_eta*z_zeta)*jacob0[idx_xi];
                xi_z0[idx_xi]=(x_eta*y_zeta-x_zeta*y_eta)*jacob0[idx_xi];
                
                eta_x0[idx_xi]=(y_zeta*z_xi-y_xi*z_zeta)*jacob0[idx_xi];
                eta_y0[idx_xi]=(x_xi*z_zeta-x_zeta*z_xi)*jacob0[idx_xi];
                eta_z0[idx_xi]=(x_zeta*y_xi-x_xi*y_zeta)*jacob0[idx_xi];
                
                zeta_x0[idx_xi]=(y_xi*z_eta-y_eta*z_xi)*jacob0[idx_xi];
                zeta_y0[idx_xi]=(x_eta*z_xi-x_xi*z_eta)*jacob0[idx_xi];
                zeta_z0[idx_xi]=(x_xi*y_eta-x_eta*y_xi)*jacob0[idx_xi];
                
                
            }
        }
    }
    
}

void nuc3d::block::interpolation_derlag(const nuc3d::Field &input,nuc3d::Field &output,int dir)
{
    int Dim_node_xi=input.getSizeX();
    int Dim_node_eta=input.getSizeY();
    int Dim_node_zeta=input.getSizeZ();
    
    int tileDim_xi=(Dim_node_xi-1)/TILE_SIZE;
    int tileDim_eta=(Dim_node_eta-1)/TILE_SIZE;
    int tileDim_zeta=(Dim_node_zeta-1)/TILE_SIZE;
    
    //std::cout<<"...";
    for(int k_tile=0;k_tile!=tileDim_zeta;k_tile++)
    {
        for(int j_tile=0;j_tile!=tileDim_eta;j_tile++)
        {
            for(int i_tile=0;i_tile!=tileDim_xi;i_tile++)
            {
                double value=0.0;
                int ibeg_center=i_tile*TILE_SIZE;
                int jbeg_center=j_tile*TILE_SIZE;
                int kbeg_center=k_tile*TILE_SIZE;
                
                int iend_center=(i_tile+1)*TILE_SIZE;
                int jend_center=(j_tile+1)*TILE_SIZE;
                int kend_center=(k_tile+1)*TILE_SIZE;
                
                int ibeg_node=i_tile*TILE_SIZE;
                int jbeg_node=j_tile*TILE_SIZE;
                int kbeg_node=k_tile*TILE_SIZE;
                
                int iend_node=(i_tile+1)*TILE_SIZE;
                int jend_node=(j_tile+1)*TILE_SIZE;
                int kend_node=(k_tile+1)*TILE_SIZE;
                
                for(int k=kbeg_center;k<kend_center;k++)
                {
                    for(int j=jbeg_center;j<jend_center;j++)
                    {
                        for(int i=ibeg_center;i<iend_center;i++)
                        {
                            value=(this->*der_lag[dir])(ibeg_node,
                                                        iend_node,
                                                        jbeg_node,
                                                        jend_node,
                                                        kbeg_node,
                                                        kend_node,
                                                        input,
                                                        i-ibeg_center,
                                                        j-jbeg_center,
                                                        k-kbeg_center);
                            output.setValue(i, j, k, value);
                            
                        }
                    }
                }
                
            }
        }
    }
}

double nuc3d::block::interpolation_derlag_center_xi(const int ibeg,
                                                    const int iend,
                                                    const int jbeg,
                                                    const int jend,
                                                    const int kbeg,
                                                    const int kend,
                                                    const Field &myfield,
                                                    int idx_xi,
                                                    int idx_eta,
                                                    int idx_zeta)
{
    double value_zeta=0.0;
    for(int k=kbeg;k<=kend;k++)
    {
        double value_eta=0.0;
        for(int j=jbeg;j<=jend;j++)
        {
            double value_xi=0.0;
            for(int i=ibeg;i<=iend;i++)
            {
                value_xi+=derlag_coeff[idx_xi][i-ibeg]*myfield.getValue(i, j, k);
            }
            value_eta+=lag_coeff[idx_eta][j-jbeg]*value_xi;
        }
        value_zeta+=lag_coeff[idx_zeta][k-kbeg]*value_eta;
    }
    
    return value_zeta;
}

double nuc3d::block::interpolation_derlag_center_eta(const int ibeg,
                                                     const int iend,
                                                     const int jbeg,
                                                     const int jend,
                                                     const int kbeg,
                                                     const int kend,
                                                     const Field &myfield,
                                                     int idx_xi,
                                                     int idx_eta,
                                                     int idx_zeta)
{
    double value_zeta=0.0;
    for(int k=kbeg;k<=kend;k++)
    {
        double value_eta=0.0;
        for(int j=jbeg;j<=jend;j++)
        {
            double value_xi=0.0;
            for(int i=ibeg;i<=iend;i++)
            {
                value_xi+=lag_coeff[idx_xi][i-ibeg]*myfield.getValue(i, j, k);
            }
            value_eta+=derlag_coeff[idx_eta][j-jbeg]*value_xi;
        }
        value_zeta+=lag_coeff[idx_zeta][k-kbeg]*value_eta;
    }
    
    return value_zeta;
}

double nuc3d::block::interpolation_derlag_center_zeta(const int ibeg,
                                                      const int iend,
                                                      const int jbeg,
                                                      const int jend,
                                                      const int kbeg,
                                                      const int kend,
                                                      const Field &myfield,
                                                      int idx_xi,
                                                      int idx_eta,
                                                      int idx_zeta)
{
    
    double value_zeta=0.0;
    for(int k=kbeg;k<=kend;k++)
    {
        double value_eta=0.0;
        for(int j=jbeg;j<=jend;j++)
        {
            double value_xi=0.0;
            for(int i=ibeg;i<=iend;i++)
            {
                value_xi+=lag_coeff[idx_xi][i-ibeg]*myfield.getValue(i, j, k);
            }
            value_eta+=lag_coeff[idx_eta][j-jbeg]*value_xi;
        }
        value_zeta+=derlag_coeff[idx_zeta][k-kbeg]*value_eta;
    }
    
    return value_zeta;
}

void nuc3d::block::initialData(int nx0,int ny0,int nz0,physicsModel &myPhy)
{
    myPDE.initPDEData3d(nx, ny, nz, myPhy.getEqNum());
    if("Euler3d"==myPhy.getMyModelName())
    {
        myFluxes=std::make_shared<EulerData3D>(nx0,ny0,nz0,myPhy.getEqNum());
    }
    else if ("EulerReactive3d"==myPhy.getMyModelName())
    {
        //        myFluxes=std::make_shared<EulerReactiveData3D>(nx0,ny0,nz0,myPhy.getEqNum());
    }
    else if ("NavierStokes3d"==myPhy.getMyModelName())
    {
        myFluxes=std::make_shared<NaiverStokesData3d>(nx0,ny0,nz0,myPhy.getEqNum());
    }
    else if ("NaiverStokesReactive3d"==myPhy.getMyModelName())
    {
        //        myFluxes=std::make_shared<NaiverStokesReactiveData3d>(nx0,ny0,nz0,myPhy.getEqNum());
    }
    else
    {
        std::cout <<"Model Name:"<<myPhy.getMyModelName()
        <<"does not exist!"
        <<std::endl;
        exit(-1);
    }
    
    for(int n=0;n<myPhy.getEqNum();n++)
    {
        mybuffer.push_back(bufferData(nx,ny,nz,bfsize));
        OutPutValue_prim.push_back(Field(nx,ny,nz));
    }
    
    for (int n=0; n<3; n++)
    {
        OutPutValue_acoust.push_back(Field(nx,ny,nz));
    }
    
    myPost=std::make_shared<postproc>(nx0,ny0,nz0);
}


void nuc3d::block::initialQ(IOController &myIO,physicsModel &myPhyMod)
{
    Field &jac=myFluxes->getJac();
    VectorField &Q=myPDE.getQ();
    
    double rho,u,v,w,p;
    double rhou,rhov,rhow,rhoe;
    double x,y,z;
    double jacobian;
    
    double gamma=myPhyMod.getGamma();
    double mach=myPhyMod.getMach();
    int fp=myIO.getStep("Benchmark");
    
    int nx0=jac.getSizeX();
    int ny0=jac.getSizeY();
    int nz0=jac.getSizeZ();
    
    for(int k=0;k<nz0;k++)
    {
        for(int j=0;j<ny0;j++)
        {
            for(int i=0;i<nx0;i++)
            {
                x=xyz_center[0].getValue(i, j, k);
                y=xyz_center[1].getValue(i, j, k);
                z=xyz_center[2].getValue(i, j, k);
                jacobian=jac.getValue(i, j, k);
                
                (this->*myInitial[fp])(rho,u,v,w,p,mach,x,y,z,gamma);
                
                rhou=rho*u;
                rhov=rho*v;
                rhow=rho*w;
                rhoe=(p/(gamma-1.0)+0.5*rho*(u*u+v*v+w*w));
                
                Q[0].setValue(i, j, k, rho/jacobian);
                Q[1].setValue(i, j, k, rhou/jacobian);
                Q[2].setValue(i, j, k, rhov/jacobian);
                Q[3].setValue(i, j, k, rhow/jacobian);
                Q[4].setValue(i, j, k, rhoe/jacobian);
            }
        }
    }
}

void nuc3d::block::solve(fieldOperator3d &myOP,
                         physicsModel &myPhyMod,
                         MPIComunicator3d_nonblocking &myMPI,
                         boundaryCondition &myBC,
                         IOController &myIO)
{
    int step=0;
    double t0=MPI_Wtime();
    while (myOP.getSteps()!=step)
    {
        myFluxes->solve(myPDE, myOP, mybuffer, myPhyMod, myMPI, myBC);
        myPDE.solve(myOP, myIO.myTimeController["cfl"],step++);
    }
    double t1=MPI_Wtime();
    
    wall_time=t1-t0;
    istep++;
    
    dt=myPDE.getDt();
    time+=dt;
    
    myIO.myTimeController["dt"]=dt;
    
    
    RES=myPDE.getRes();
}

void nuc3d::block::printStatus()
{
    std::cout<<std::setprecision(6)<<"========step = "<<istep<< "\n time = "<<time<<", dt = "<<dt<<", residual = "<<RES<<", CPU time = "<<wall_time<<"(s) "<<std::endl;
}

void nuc3d::block::Post(fieldOperator3d &myOp,
                        physicsModel &myPhys,
                        MPIComunicator3d_nonblocking &myMPI,
                        boundaryCondition &myBC,
                        IOController &myIO)
{
    Field &jac=myFluxes->getJac();
    myPhys.getPrim(jac, myPDE.getQ(), OutPutValue_prim, OutPutValue_acoust);
    
    myPost->solvePost(OutPutValue_prim,
                      OutPutValue_acoust,
                      xyz,
                      myFluxes->getXi_xyz(),
                      myFluxes->getEta_xyz(),
                      myFluxes->getZeta_xyz(),
                      myOp, mybuffer, myMPI, myBC,myIO,istep,time);
    
    if(("yes")==(myIO.getType("PostProc")))
    {
        myPost->OutputPost(OutPutValue_prim,
                           OutPutValue_acoust,
                           xyz,
                           myFluxes->getXi_xyz(),
                           myFluxes->getEta_xyz(),
                           myFluxes->getZeta_xyz(),
                           myOp, mybuffer, myMPI, myBC,myIO,istep,time);
        
    }
    
    
}

void nuc3d::block::Output(fieldOperator3d &myOp,
                          physicsModel &myPhys,
                          MPIComunicator3d_nonblocking &myMPI,
                          boundaryCondition &myBC,
                          IOController &myIO)
{
    Field &jac=myFluxes->getJac();
    myPhys.getPrim(jac, myPDE.getQ(), OutPutValue_prim, OutPutValue_acoust);
    
    if(("yes")==(myIO.getType("Binary")))
        outputQ_binary(myMPI.getMyId(),myPhys);
    
    if(("yes")==(myIO.getType("Tecplot")))
    {
        outputQ_tecplot(myMPI.getMyId(),myPhys);
        
    }
    
}


void nuc3d::block::initial_default(double &rho,double &u,double &v,double &w,double &p,double &mach,double &x,double &y,double &z,double &gamma)
{
    rho=1.0;
    u=1.0;
    v=0.0;
    w=0.0;
    p=1.0/(mach*mach*gamma);
    
}

void nuc3d::block::initial_ivc(double &rho,double &u,double &v,double &w,double &p,double &mach,double &x,double &y,double &z,double &gamma)
{
    double pie=4.0*atan(1.0);
    double b=0.5;
    double x_c=5.0;
    double y_c=5.0;
    
    double r=std::sqrt(std::pow(x-x_c,2)+std::pow(y-y_c,2));
    rho=std::pow((1-(gamma-1.0)*b*b*std::exp(1-r*r)/(8.0*gamma*pie*pie)),2.5);
    u=0.5-b/(2.0*pie)*exp((1-r*r)/2)*(y-y_c);
    v=b/(2.0*pie)*exp((1-r*r)/2)*(x-x_c);
    w=0.0;
    p=std::pow(rho,gamma);
    
}
void nuc3d::block::initial_taylorgreen(double &rho,double &u,double &v,double &w,double &p,double &mach,double &x,double &y,double &z,double &gamma)
{
    rho=1.0;
    u=std::sin(x)*std::cos(y)*std::cos(z);
    v=-std::cos(x)*std::sin(y)*std::cos(z);
    w=0.0;
    p=1.0/(gamma*mach*mach)+((std::cos(2.0*z)+2.0)*(std::cos(2.0*x)+std::cos(2.0*y))-2.0)/16.0;
}

void nuc3d::block::outputQ_tecplot(int myID,physicsModel &myPhys)
{
    std::string forename_flow = ("flowData/NUC3d_step_");
    std::string step;
    std::string mid("_id_");
    std::string id;
    std::string tailname = (".dat");
    
    std::stringstream ss_step,ss_id;
    ss_step << istep;
    ss_step >> step;
    
    ss_id<<myID;
    ss_id>>id;
    
    std::string filename_flow = forename_flow + step + mid + id + tailname;
    
    std::ofstream myIOfile;
    
    myIOfile.open(filename_flow);
    
    std::string TECplotHeader[2]={"title=NUC3d\n",
        "variables=x,y,z"};
    
    myIOfile<<TECplotHeader[0]
    <<TECplotHeader[1];
    
    for(int i=0;i<(OutPutValue_prim.size()+OutPutValue_acoust.size());i++)
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
    <<xyz.size()+OutPutValue_prim.size()+OutPutValue_acoust.size()
    <<"]=CELLCENTERED)\n";
    
    for(auto iter=xyz.begin();iter!=xyz.end();iter++)
    {
        writeField(myIOfile, *iter);
    }
    
    for(auto iter=OutPutValue_prim.begin();iter!=OutPutValue_prim.end();iter++)
    {
        writeField(myIOfile, *iter);
    }
    
    for(auto iter=OutPutValue_acoust.begin();iter!=OutPutValue_acoust.end();iter++)
    {
        writeField(myIOfile, *iter);
    }
    
    myIOfile.close();
    
}

void nuc3d::block::outputQ_binary(int myID,physicsModel &myPhys)
{
    std::string forename_flow = ("flowData/NUC3d_step_");
    std::string step;
    std::string mid("_ID_");
    std::string id;
    std::string tailname = (".bin");
    
    std::stringstream ss_step,ss_id;
    ss_step << istep;
    ss_step >> step;
    
    ss_id<<myID;
    ss_id>>id;
    
    std::string filename_flow = forename_flow + step + mid + id + tailname;
    
    std::ofstream myIOfile;
    
    myIOfile.open(filename_flow,std::ios::out|std::ios::binary);
    VectorField &Q=myPDE.getQ();
    myIOfile.write(reinterpret_cast<char *>(&time), sizeof(time));
    for(auto iter=Q.begin();iter!=Q.end();iter++)
    {
        writeField_binary(myIOfile,*iter);
    }
    
    myIOfile.close();
    
}

void nuc3d::block::inputQ_binary(int myID,int step0)
{
    std::string forename_flow = ("flowData/NUC3d_step_");
    std::string step;
    std::string mid("_ID_");
    std::string id;
    std::string tailname = (".bin");
    
    std::stringstream ss_step,ss_id;
    ss_step << step0;
    ss_step >> step;
    
    ss_id<<myID;
    ss_id>>id;
    
    std::string filename_flow = forename_flow + step + mid + id + tailname;
    
    std::ifstream myIOfile;
    
    myIOfile.open(filename_flow,std::ios::in|std::ios::binary);
    VectorField &Q=myPDE.getQ();
    
    myIOfile.read(reinterpret_cast<char *>(&time), sizeof (time));
    istep=step0;
    for(auto iter=Q.begin();iter!=Q.end();iter++)
    {
        readField_binary(myIOfile,*iter);
    }
    
    myIOfile.close();
    
}

void nuc3d::block::writeField_binary(std::ofstream &myFile, nuc3d::Field &myField)
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
                myFile.write(reinterpret_cast<char *>(&value), sizeof(value));
            }
        }
    }
}


void nuc3d::block::readField_binary(std::ifstream &myFile, nuc3d::Field &myField)
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
                double value;
                myFile.read(reinterpret_cast<char *>(&value), sizeof(value));
                myField.setValue(i,j,k,value);
            }
        }
    }
}

void nuc3d::block::outputGEO_tecplot(int myID)
{
    std::string forename_flow = ("flowData/GEO_");
    std::string mid("_id_");
    std::string id;
    std::string tailname = (".dat");
    
    std::stringstream ss_id;
    
    ss_id<<myID;
    ss_id>>id;
    
    std::string filename_flow = forename_flow + mid + id + tailname;
    
    std::ofstream myIOfile;
    
    myIOfile.open(filename_flow);
    
    std::string TECplotHeader[2]={"title=NUC3d\n",
        "variables=x,y,z"};
    
    myIOfile<<TECplotHeader[0]
    <<TECplotHeader[1];
    
    for(int i=0;i<(myFluxes->getXi_xyz().size()
                   +myFluxes->getEta_xyz().size()
                   +myFluxes->getZeta_xyz().size());i++)
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
    <<xyz.size()+myFluxes->getXi_xyz().size()
    +myFluxes->getEta_xyz().size()
    +myFluxes->getZeta_xyz().size()
    <<"]=CELLCENTERED)\n";
    
    for(auto iter=xyz.begin();iter!=xyz.end();iter++)
    {
        writeField(myIOfile, *iter);
    }
    
    for(auto iter=myFluxes->getXi_xyz().begin();iter!=myFluxes->getXi_xyz().end();iter++)
    {
        writeField(myIOfile, *iter);
    }
    
    for(auto iter=myFluxes->getEta_xyz().begin();iter!=myFluxes->getEta_xyz().end();iter++)
    {
        writeField(myIOfile, *iter);
    }
    
    for(auto iter=myFluxes->getZeta_xyz().begin();iter!=myFluxes->getZeta_xyz().end();iter++)
    {
        writeField(myIOfile, *iter);
    }
    myIOfile.close();
    
}

void nuc3d::block::solve_grad(Field &myField,
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
    
    SolveGradEta(feta,myOP,myBf,myMPI,deta,fdID,typeL,typeR);
    
    typeL=myBC.getBCtype(4);
    typeR=myBC.getBCtype(5);
    
    solveGradZeta(fzeta,myOP,myBf,myMPI,dzeta,fdID,typeL,typeR);
}

void nuc3d::block::solveGradXi(Field &myField,
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

void nuc3d::block::SolveGradEta(Field &myField,
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

void nuc3d::block::solveGradZeta(Field &myField,
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

