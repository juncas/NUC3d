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
    
    std::string forename_mesh = ("mesh_");
    std::string forename_flow = ("flow_");
    std::string midname;
    std::string tailname = (".x");
    
    ss << ProcId;
    ss >> midname;
    
    std::string filename_mesh = forename_mesh + midname + tailname;
    std::string filename_flow = forename_flow + midname + tailname;
    
    std::ifstream myFile;
    std::ofstream myFile_o("IOtest.dat");
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
        std::cout<<"Start reading mesh data..."<<std::endl;
        for(auto iter=xyz.begin();iter!=xyz.end();iter++)
            readField(myFile,*iter);
        std::cout<<"Mesh data has been read!"<<std::endl;
        
        std::cout<<"Start calculating mesh data..."<<std::endl;
        getXYZ_center();
        std::cout<<"Center data has been calculated..."<<std::endl;
        std::cout<<"Mesh data has been recalculated!"<<std::endl;
        //writeField(myFile_o,*iter);
    }
    else
    {
        std::cout<<"File \'"<<filename_mesh<<"\' does not exist!"<<std::endl;
        exit(-1);
    }
    myFile.close();
    
    initialData(nx, ny, nz, myPhyMod);
    getJacobians();
    
//    for(auto iter=myFluxes->getXi_xyz().begin();iter!=myFluxes->getXi_xyz().end();iter++)
//        writeField(myFile_o, *iter);
//    for(auto iter=myFluxes->getEta_xyz().begin();iter!=myFluxes->getEta_xyz().end();iter++)
//            writeField(myFile_o, *iter);
    for(auto iter=myFluxes->getZeta_xyz().begin();iter!=myFluxes->getZeta_xyz().end();iter++)
        writeField(myFile_o, *iter);
    
    std::cout<<"Jacobians has been calculated..."<<std::endl;
    initialQ();
    
    //for(auto iter=myPDE.getQ().begin();iter!=myPDE.getQ().end();iter++)
    //    writeField(myFile_o, *iter);
    outputQ();
    myBC.initialBC(mybuffer,myMPI);
    
}
void nuc3d::block::writeField(std::ofstream &myFile, nuc3d::Field &myField)
{
    int nx0=myField.getSizeX();
    int ny0=myField.getSizeY();
    int nz0=myField.getSizeZ();
    myFile<<nx0<<" "<<ny0<<" "<<nz0<<"\n";
    for(int k=0;k<nz0;k++)
    {
        for(int j=0;j<ny0;j++)
        {
            for(int i=0;i<nx0;i++)
            {
                double value=myField.getValue(i,j,k);
                myFile<<std::setprecision(12)<<value<<" ";
            }
        }
    }
    myFile<<"\n";
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

void nuc3d::block::getJacobians()
{
    VectorField &xi_xyz=myFluxes->getXi_xyz(); // xyz_Xi
    VectorField &eta_xyz=myFluxes->getEta_xyz();// xyz_Eta
    VectorField &zeta_xyz=myFluxes->getZeta_xyz();// xyz_Zeta
    
    Field &jac=myFluxes->getJac();
    
    std::cout<<"Start calculate xyz_xi, xyz_eta, xyz_zta ..."<<std::endl;
    
    
    auto beg=xyz.begin();
    auto end=xyz.end();
    
    for(auto iter=beg;iter!=end;iter++)
    {
        interpolation_derlag(*iter, xi_xyz[iter-beg],0);
        interpolation_derlag(*iter, eta_xyz[iter-beg],1);
        interpolation_derlag(*iter, zeta_xyz[iter-beg],2);
    }
    
    std::cout<<std::endl;
    
    std::cout<<"xx_xyz has been calculated ..."<<std::endl;
    std::cout<<"Calculating Jacobians calculated ..."<<std::endl;
    
    
    int nx0=jac.getSizeX();
    int ny0=jac.getSizeY();
    int nz0=jac.getSizeZ();
    for(int k=0;k<nz0;k++)
    {
        for(int j=0;j<ny0;j++)
        {
            for(int i=0;i<nx0;i++)
            {
                double x_xi=xi_xyz[0].getValue(i, j, k);
                double y_xi=xi_xyz[1].getValue(i, j, k);
                double z_xi=xi_xyz[2].getValue(i, j, k);
                
                double x_eta=eta_xyz[0].getValue(i, j, k);
                double y_eta=eta_xyz[1].getValue(i, j, k);
                double z_eta=eta_xyz[2].getValue(i, j, k);
                
                double x_zeta=zeta_xyz[0].getValue(i, j, k);
                double y_zeta=zeta_xyz[1].getValue(i, j, k);
                double z_zeta=zeta_xyz[2].getValue(i, j, k);
                
                double jacob=1.0/(x_xi*(y_eta*z_zeta-y_zeta*z_eta)
                                  +x_eta*(y_zeta*z_xi-y_xi*z_zeta)
                                  +x_zeta*(y_xi*z_eta-y_eta*z_xi));
                
                double xi_x=(y_eta*z_zeta-y_zeta*z_eta)*jacob;
                double xi_y=(y_zeta*z_xi-y_xi*z_zeta)*jacob;
                double xi_z=(y_xi*z_eta-y_eta*z_xi)*jacob;
                
                double eta_x=(x_zeta*z_eta-x_eta*z_zeta)*jacob;
                double eta_y=(x_xi*z_zeta-x_zeta*z_xi)*jacob;
                double eta_z=(x_eta*z_xi-x_xi*z_eta)*jacob;
                
                double zeta_x=(x_eta*y_zeta-x_zeta*y_eta)*jacob;
                double zeta_y=(x_zeta*y_xi-x_xi*y_zeta)*jacob;
                double zeta_z=(x_xi*y_eta-x_eta*y_xi)*jacob;
                
                xi_xyz[0].setValue(i, j, k, xi_x);
                xi_xyz[1].setValue(i, j, k, xi_y);
                xi_xyz[2].setValue(i, j, k, xi_z);
                
                eta_xyz[0].setValue(i, j, k, eta_x);
                eta_xyz[1].setValue(i, j, k, eta_y);
                eta_xyz[2].setValue(i, j, k, eta_z);
                
                zeta_xyz[0].setValue(i, j, k, zeta_x);
                zeta_xyz[1].setValue(i, j, k, zeta_y);
                zeta_xyz[2].setValue(i, j, k, zeta_z);
                
                jac.setValue(i, j, k, jacob);
                
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
    
    std::cout<<"...";
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
    else if ("NaiverStokes3d"==myPhy.getMyModelName())
    {
        //        myFluxes=std::make_shared<NaiverStokesData3d>(nx0,ny0,nz0,myPhy.getEqNum());
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
        mybuffer.push_back(bufferData(nx,ny,nz,bfsize));
    
    std::cout<<"Field data has been allocated!"<<std::endl;
}

std::string TECplotHeader("title=NUC3d\nvariables=\"x\",\"y\",\"z\",\"rho\",\"u\",\"v\",\"w\",\"p\"\n");

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
        myPDE.solve(myOP, step++);
    }
    double t1=MPI_Wtime();
    
    wall_time=t1-t0;
    istep++;
    myIO.myIOController["currentStep"]=istep;
    
    dt=myPDE.getDt();
    time+=dt;
    
    myIO.myTimeController["dt"]=dt;
    myIO.myTimeController["currentTime"]=time;
    
    
    RES=myPDE.getRes();
}

void nuc3d::block::printStatus()
{
    std::cout<<std::setprecision(6)<<"========step = "<<istep<< "\n time = "<<time<<", dt = "<<dt<<", residual = "<<RES<<", CPU time = "<<wall_time<<"(s) "<<std::endl;
}

void nuc3d::block::initialQ()
{
    Field &jac=myFluxes->getJac();
    VectorField &Q=myPDE.getQ();
    
    double rho,u,v,w,p;
    double rhou,rhov,rhow,rhoe;
    double x,y,z;
    double jacobian;
    
    double x_c=5.0;
    double y_c=5.0;
    double gamma=1.4;
    double b=0.5;
    double pie=4.0*atan(1.0);
    
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
                
                double r_sq=std::pow((x-x_c),2)+std::pow((y-y_c), 2);
                
                rho=pow(1.0-(gamma-1.0)*b*b/8.0/gamma/pie/pie*std::exp(1.0-r_sq),1/(gamma-1.0));
                
                u=0.5-b/2.0/pie*exp((1-r_sq)/2.0)*(y-y_c);
                v=b/2.0/pie*exp((1-r_sq)/2.0)*(x-x_c);
                w=0.0;
                
                p=pow(rho,gamma);
                
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


void nuc3d::block::outputQ()
{
    Field &jac=myFluxes->getJac();
    VectorField &Q=myPDE.getQ();
    
    double rho,u,v,w,p;
    double x,y,z;
    double jacobian;
    double gamma=1.4;
    
    int nx0=jac.getSizeX();
    int ny0=jac.getSizeY();
    int nz0=jac.getSizeZ();
    
    std::string forename_flow = ("flow_");
    std::string midname;
    std::string tailname = (".dat");
    
    std::stringstream ss;
    ss << istep;
    ss >> midname;
    
    std::string filename_flow = forename_flow + midname + tailname;
    

    std::ofstream myIOfile;
    
    myIOfile.open(filename_flow);
    myIOfile<<TECplotHeader;
    myIOfile<<"Zone I = "<<nx0<<", J= "<<ny0<<", K="<<nz0<<"\n";
    
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
                
                rho=Q[0].getValue(i, j, k)*jacobian;
                u=Q[1].getValue(i, j, k)*jacobian;
                v=Q[2].getValue(i, j, k)*jacobian;
                w=Q[3].getValue(i, j, k)*jacobian;
                p=(Q[4].getValue(i, j, k)*jacobian-0.5*rho*(u*u+v*v+w*w))*(gamma-1.0);
                
                myIOfile<<x<<" "<<y<<" "<<z<<" "<<rho<<" "<<u<<" "<<v<<" "<<w<<" "<<p<<"\n";
               
            }
        }
    }
    
}

