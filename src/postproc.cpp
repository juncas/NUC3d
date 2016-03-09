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
#include "physicsModel.h"
#include "fieldOperator.h"

nuc3d::postproc::postproc(int nx0,int ny0,int nz0):
nx(nx0),
ny(ny0),
nz(nz0),
temp_xi(nx0,ny0,nz0),
temp_eta(ny0,nz0,nx0),
temp_zeta(nz0,nx0,ny0),
dtempdxi(nx0,ny0,nz0),
dtempdeta(ny0,nz0,nx0),
dtempdzeta(nz0,nx0,ny0)
{
    
    std::string filename("inp/PostProc.in");
    std::ifstream file(filename);
    
    std::string word0;
    std::string word1;
    if(file)
    {
        while(file>>word0>>word1)
        {
            std::istringstream value0(word1);
            int value;
            value0>>value;
            
            variableScalarInt.push_back(value);
        }
    }
    else
    {
        std::cout<<"File inp/PostProc.in does not exist!"<<std::endl;
        exit(-1);
    }
    
    for(int n=0;n<variableScalarInt[0];n++)
    {
        variableField.push_back(Field(nx,ny,nz,0.0));
    }
    
    for(int n=0;n<variableScalarInt[1];n++)
    {
        variableScalarDouble.push_back(0.0);
    }
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
                                physicsModel myPhys,
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
    
    solveGrad(u, myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(variableField[0],variableField[1],variableField[2],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    solveGrad(v, myOP, myBf[1], myMPI,myBC,1,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(variableField[3],variableField[4],variableField[5],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    solveGrad(w, myOP, myBf[2], myMPI,myBC,2,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(variableField[6],variableField[7],variableField[8],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
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

void nuc3d::postproc::solveTemporal(VectorField &prims,
                                    VectorField &acous,
                                    VectorField &xi_xyz,
                                    VectorField &eta_xyz,
                                    VectorField &zeta_xyz,
                                    physicsModel myPhys,
                                    fieldOperator3d &myOP,
                                    VectorBuffer &myBf,
                                    MPIComunicator3d_nonblocking &myMPI,
                                    boundaryCondition &myBC)
{
    Field &rho=prims[0];
    Field &u=prims[1];
    Field &v=prims[2];
    Field &w=prims[3];
    Field &p=prims[4];
    Field &T=acous[0];
    
    solveGrad(u, myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(variableField[0],variableField[1],variableField[2],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    solveGrad(v, myOP, myBf[1], myMPI,myBC,1,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(variableField[3],variableField[4],variableField[5],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    solveGrad(w, myOP, myBf[2], myMPI,myBC,2,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(variableField[6],variableField[7],variableField[8],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    solveGrad(p, myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(variableField[9],variableField[10],variableField[11],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    myPhys.getMiu(T, variableField[12], variableField[13]);
    
    double *prho=rho.getDataPtr();
    double *pu=u.getDataPtr();
    double *pv=v.getDataPtr();
    double *pw=w.getDataPtr();
    double *pp=p.getDataPtr();
    double *pT=T.getDataPtr();
    
    double* ux=variableField[0].getDataPtr();
    double* uy=variableField[1].getDataPtr();
    double* uz=variableField[2].getDataPtr();
    double* vx=variableField[3].getDataPtr();
    double* vy=variableField[4].getDataPtr();
    double* vz=variableField[5].getDataPtr();
    double* wx=variableField[6].getDataPtr();
    double* wy=variableField[7].getDataPtr();
    double* wz=variableField[8].getDataPtr();
    double* pmiu=variableField[12].getDataPtr();
    
    double *pOmega0=variableField[14].getDataPtr();
    double *pOmega1=variableField[15].getDataPtr();
    double *pOmega2=variableField[16].getDataPtr();
    double *pQ=variableField[17].getDataPtr();
    
    double *ptau_xx=variableField[18].getDataPtr();
    double *ptau_xy=variableField[19].getDataPtr();
    double *ptau_xz=variableField[20].getDataPtr();
    double *ptau_yy=variableField[21].getDataPtr();
    double *ptau_yz=variableField[22].getDataPtr();
    double *ptau_zz=variableField[23].getDataPtr();
    
    double *prhou=variableField[24].getDataPtr();
    double *prhov=variableField[25].getDataPtr();
    double *prhow=variableField[26].getDataPtr();
    
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                
                pOmega0[idx_xi]=wy[idx_xi]-vz[idx_xi];
                pOmega1[idx_xi]=uz[idx_xi]-wx[idx_xi];
                pOmega2[idx_xi]=vx[idx_xi]-uy[idx_xi];
                
                pQ[idx_xi]=ux[idx_xi]*vy[idx_xi]
                +ux[idx_xi]*wz[idx_xi]
                +vy[idx_xi]*wz[idx_xi]
                -uy[idx_xi]*vx[idx_xi]
                -uz[idx_xi]*wx[idx_xi]
                -vz[idx_xi]*wy[idx_xi];
                
                double VelGrad=ux[idx_xi]+vy[idx_xi]+wz[idx_xi];
                
                double miu0=pMiu[idx_xi];
                
                ptau_xx[idx_xi]=miu0*(2.0*ux[idx_xi]-2.0/3.0*VelGrad);
                ptau_xy[idx_xi]=miu0*(uy[idx_xi]+vx[idx_xi]-2.0/3.0*VelGrad);
                ptau_xz[idx_xi]=miu0*(uz[idx_xi]+wx[idx_xi]-2.0/3.0*VelGrad);
                ptau_yy[idx_xi]=miu0*(2.0*vy[idx_xi]-2.0/3.0*VelGrad);
                ptau_yz[idx_xi]=miu0*(vz[idx_xi]+wy[idx_xi]-2.0/3.0*VelGrad);
                ptau_zz[idx_xi]=miu0*(2.0*wz[idx_xi]-2.0/3.0*VelGrad);
                
                prhou=pu[idx_xi]*prho[idx_xi];
                prhov=pv[idx_xi]*prho[idx_xi];
                prhow=pw[idx_xi]*prho[idx_xi];
            }
        }
    }
    
    solveGrad(variableField[18], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(variableField[27],variableField[28],variableField[29],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    solveGrad(variableField[19], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(variableField[30],variableField[31],variableField[32],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    solveGrad(variableField[20], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(variableField[33],variableField[34],variableField[35],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    solveGrad(variableField[21], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(variableField[36],variableField[37],variableField[38],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    solveGrad(variableField[22], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(variableField[39],variableField[40],variableField[41],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    solveGrad(variableField[23], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(variableField[42],variableField[43],variableField[44],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    solveGrad(variableField[24], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(variableField[45],variableField[46],variableField[47],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    solveGrad(variableField[25], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(variableField[48],variableField[49],variableField[50],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    solveGrad(variableField[26], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(variableField[51],variableField[52],variableField[53],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    
}

void nuc3d::postproc::solveAveraged(VectorField &prims,
                                    VectorField &acous,
                                    VectorField &xi_xyz,
                                    VectorField &eta_xyz,
                                    VectorField &zeta_xyz,
                                    fieldOperator3d &myOP,
                                    VectorBuffer &myBf,
                                    MPIComunicator3d_nonblocking &myMPI,
                                    boundaryCondition &myBC)
{
    
    double *prho=rho.getDataPtr();
    double *pu=u.getDataPtr();
    double *pv=v.getDataPtr();
    double *pw=w.getDataPtr();
    double *pp=p.getDataPtr();
    double *pT=T.getDataPtr();
    
    double* ux=variableField[0].getDataPtr();
    double* uy=variableField[1].getDataPtr();
    double* uz=variableField[2].getDataPtr();
    double* vx=variableField[3].getDataPtr();
    double* vy=variableField[4].getDataPtr();
    double* vz=variableField[5].getDataPtr();
    double* wx=variableField[6].getDataPtr();
    double* wy=variableField[7].getDataPtr();
    double* wz=variableField[8].getDataPtr();
    double* ppx=variableField[9].getDataPtr();
    double* ppy=variableField[10].getDataPtr();
    double* ppz=variableField[11].getDataPtr();
    double* pmiu=variableField[12].getDataPtr();
    
    double *pOmega0=variableField[14].getDataPtr();
    double *pOmega1=variableField[15].getDataPtr();
    double *pOmega2=variableField[16].getDataPtr();
    double *pQ=variableField[17].getDataPtr();
    
    double *ptau_xx=variableField[18].getDataPtr();
    double *ptau_xy=variableField[19].getDataPtr();
    double *ptau_xz=variableField[20].getDataPtr();
    double *ptau_yy=variableField[21].getDataPtr();
    double *ptau_yz=variableField[22].getDataPtr();
    double *ptau_zz=variableField[23].getDataPtr();
    
    double *prhou=variableField[24].getDataPtr();
    double *prhov=variableField[25].getDataPtr();
    double *prhow=variableField[26].getDataPtr();

    
    double *pOmega0=variableField[14].getDataPtr();
    double *pOmega1=variableField[15].getDataPtr();
    double *pOmega2=variableField[16].getDataPtr();
    double *pQ=variableField[17].getDataPtr();
    
    double *ptau_xx=variableField[18].getDataPtr();
    double *ptau_xy=variableField[19].getDataPtr();
    double *ptau_xz=variableField[20].getDataPtr();
    double *ptau_yy=variableField[21].getDataPtr();
    double *ptau_yz=variableField[22].getDataPtr();
    double *ptau_zz=variableField[23].getDataPtr();
    
    double *prhou=variableField[24].getDataPtr();
    double *prhov=variableField[25].getDataPtr();
    double *prhow=variableField[26].getDataPtr();
    
    double *pdrhoudx=variableField[45].getDataPtr();
    double *pdrhoudy=variableField[46].getDataPtr();
    double *pdrhoudz=variableField[47].getDataPtr();
    
    double *pdrhovdx=variableField[48].getDataPtr();
    double *pdrhovdy=variableField[49].getDataPtr();
    double *pdrhovdz=variableField[50].getDataPtr();
    
    double *pdrhowdx=variableField[51].getDataPtr();
    double *pdrhowdy=variableField[52].getDataPtr();
    double *pdrhowdz=variableField[53].getDataPtr();
    
    double *paver_rho=variableField[54].getDataPtr();
    double *paver_u=variableField[55].getDataPtr();
    double *paver_v=variableField[56].getDataPtr();
    double *paver_w=variableField[57].getDataPtr();
    double *paver_p=variableField[58].getDataPtr();
    double *paver_T=variableField[59].getDataPtr();
    
    //    //TKE budaget terms
    //    double *pTKE=variableField[60].getDataPtr();
    //
    //    //
    //    //Production
    //    double *pP_tke=variableField[61].getDataPtr();
    
    double *paver_rhou=variableField[62].getDataPtr();
    double *paver_rhov=variableField[63].getDataPtr();
    double *paver_rhow=variableField[64].getDataPtr();
    
    double *paver_rhouu=variableField[65].getDataPtr();
    double *paver_rhouv=variableField[66].getDataPtr();
    double *paver_rhouw=variableField[67].getDataPtr();
    double *paver_rhovv=variableField[68].getDataPtr();
    double *paver_rhovw=variableField[69].getDataPtr();
    double *paver_rhoww=variableField[70].getDataPtr();
    
    double *pdaver_rhoudx=variableField[71].getDataPtr();
    double *pdaver_rhoudy=variableField[72].getDataPtr();
    double *pdaver_rhoudz=variableField[73].getDataPtr();
    
    double *pdaver_rhovdx=variableField[74].getDataPtr();
    double *pdaver_rhovdy=variableField[75].getDataPtr();
    double *pdaver_rhovdz=variableField[76].getDataPtr();
    
    double *pdaver_rhowdx=variableField[77].getDataPtr();
    double *pdaver_rhowdy=variableField[78].getDataPtr();
    double *pdaver_rhowdz=variableField[79].getDataPtr();
    
    //    double *pT_tke=variableField[80].getDataPtr();
    //    double *paver_Tx=variableField[81].getDataPtr();
    //    double *paver_Ty=variableField[82].getDataPtr();
    //    double *paver_Tz=variableField[83].getDataPtr();
    
    double *paver_rhouuu=variableField[84].getDataPtr();
    double *paver_rhouuv=variableField[85].getDataPtr();
    double *paver_rhouuw=variableField[86].getDataPtr();
    double *paver_rhovvu=variableField[87].getDataPtr();
    double *paver_rhovvv=variableField[88].getDataPtr();
    double *paver_rhovvw=variableField[89].getDataPtr();
    double *paver_rhowwu=variableField[90].getDataPtr();
    double *paver_rhowwv=variableField[91].getDataPtr();
    double *paver_rhowww=variableField[92].getDataPtr();
    //    double *pdTdx=variableField[93].getDataPtr();
    //    double *pdTdy=variableField[94].getDataPtr();
    //    double *pdTdz=variableField[95].getDataPtr();
    
    double *pPIE_tke=variableField[96].getDataPtr();
    double *pPIE_tke_0=variableField[97].getDataPtr();
    double *pPIE_tke_1=variableField[98].getDataPtr();
    double *pPIE_tke_2=variableField[99].getDataPtr();
    
    double *paver_dpdx=variableField[100].getDataPtr();
    double *paver_dpdy=variableField[101].getDataPtr();
    double *paver_dpdz=variableField[102].getDataPtr();
    
    double *paver_pu=variableField[103].getDataPtr();
    double *paver_pv=variableField[104].getDataPtr();
    double *paver_pw=variableField[105].getDataPtr();
    
    double *paver_pdudx=variableField[106].getDataPtr();
    double *paver_pdvdy=variableField[107].getDataPtr();
    double *paver_pdwdz=variableField[108].getDataPtr();
    
    double *paver_dudx=variableField[109].getDataPtr();
    double *paver_dvdy=variableField[110].getDataPtr();
    double *paver_dwdz=variableField[111].getDataPtr();
    
    double *paver_dpudx=variableField[112].getDataPtr();
    double *paver_dpvdy=variableField[113].getDataPtr();
    double *paver_dpwdz=variableField[114].getDataPtr();
    
    double *paver_pfuf=variableField[115].getDataPtr();
    double *paver_pfvf=variableField[116].getDataPtr();
    double *paver_pfwf=variableField[117].getDataPtr();
    
    double *pM_tke=variableField[118].getDataPtr();
    
    double *paver_tau_xx=variableField[119].getDataPtr();
    double *paver_tau_xy=variableField[120].getDataPtr();
    double *paver_tau_xz=variableField[121].getDataPtr();
    double *paver_tau_yy=variableField[122].getDataPtr();
    double *paver_tau_yz=variableField[123].getDataPtr();
    double *paver_tau_zz=variableField[124].getDataPtr();
    //
    double *paver_dtau_xxdx=variableField[125].getDataPtr();
    double *paver_dtau_xydy=variableField[126].getDataPtr();
    double *paver_dtau_xzdz=variableField[127].getDataPtr();
    
    double *paver_dtau_yxdx=variableField[128].getDataPtr();
    double *paver_dtau_yydy=variableField[129].getDataPtr();
    double *paver_dtau_yzdz=variableField[130].getDataPtr();
    
    double *paver_dtau_zxdx=variableField[131].getDataPtr();
    double *paver_dtau_zydy=variableField[132].getDataPtr();
    double *paver_dtau_zzdz=variableField[133].getDataPtr();


    //TKE and Production
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                
                paver_rho[idx_xi]=(paver_rho[idx_xi]*averStep+prho[idx_xi])/(averStep+1.0);
                paver_u[idx_xi]=(paver_u[idx_xi]*averStep+pu[idx_xi])/(averStep+1.0);
                paver_v[idx_xi]=(paver_v[idx_xi]*averStep+pv[idx_xi])/(averStep+1.0);
                paver_w[idx_xi]=(paver_w[idx_xi]*averStep+pw[idx_xi])/(averStep+1.0);
                paver_p[idx_xi]=(paver_p[idx_xi]*averStep+pp[idx_xi])/(averStep+1.0);
                paver_T[idx_xi]=(paver_T[idx_xi]*averStep+pT[idx_xi])/(averStep+1.0);
                
                //
                
                paver_rhou[idx_xi]=(paver_rhou[idx_xi]*averStep+prho[idx_xi]*pu[idx_xi])/(averStep+1.0);
                paver_rhov[idx_xi]=(paver_rhov[idx_xi]*averStep+prho[idx_xi]*pv[idx_xi])/(averStep+1.0);
                paver_rhow[idx_xi]=(paver_rhow[idx_xi]*averStep+prho[idx_xi]*pw[idx_xi])/(averStep+1.0);
                
                paver_rhouu[idx_xi]=(paver_rhouu[idx_xi]*averStep+prho[idx_xi]*pu[idx_xi]*pu[idx_xi])/(averStep+1.0);
                paver_rhouv[idx_xi]=(paver_rhouv[idx_xi]*averStep+prho[idx_xi]*pu[idx_xi]*pv[idx_xi])/(averStep+1.0);
                paver_rhouw[idx_xi]=(paver_rhouw[idx_xi]*averStep+prho[idx_xi]*pu[idx_xi]*pw[idx_xi])/(averStep+1.0);
                paver_rhovv[idx_xi]=(paver_rhovv[idx_xi]*averStep+prho[idx_xi]*pv[idx_xi]*pv[idx_xi])/(averStep+1.0);
                paver_rhovw[idx_xi]=(paver_rhovw[idx_xi]*averStep+prho[idx_xi]*pv[idx_xi]*pw[idx_xi])/(averStep+1.0);
                paver_rhoww[idx_xi]=(paver_rhoww[idx_xi]*averStep+prho[idx_xi]*pw[idx_xi]*pw[idx_xi])/(averStep+1.0);
                
                pdaver_rhoudx[idx_xi]=(pdaver_rhoudx[idx_xi]*averStep+pdrhoudx[idx_xi])/(averStep+1.0);
                pdaver_rhoudy[idx_xi]=(pdaver_rhoudy[idx_xi]*averStep+pdrhoudy[idx_xi])/(averStep+1.0);
                pdaver_rhoudz[idx_xi]=(pdaver_rhoudz[idx_xi]*averStep+pdrhoudz[idx_xi])/(averStep+1.0);
                
                pdaver_rhovdx[idx_xi]=(pdaver_rhovdx[idx_xi]*averStep+pdrhovdx[idx_xi])/(averStep+1.0);
                pdaver_rhovdy[idx_xi]=(pdaver_rhovdy[idx_xi]*averStep+pdrhovdy[idx_xi])/(averStep+1.0);
                pdaver_rhovdz[idx_xi]=(pdaver_rhovdz[idx_xi]*averStep+pdrhovdz[idx_xi])/(averStep+1.0);
                
                pdaver_rhowdx[idx_xi]=(pdaver_rhowdx[idx_xi]*averStep+pdrhowdx[idx_xi])/(averStep+1.0);
                pdaver_rhowdy[idx_xi]=(pdaver_rhowdy[idx_xi]*averStep+pdrhowdy[idx_xi])/(averStep+1.0);
                pdaver_rhowdz[idx_xi]=(pdaver_rhowdz[idx_xi]*averStep+pdrhowdz[idx_xi])/(averStep+1.0);
                
                paver_rhouuu[idx_xi]=(paver_rhouuu[idx_xi]*averStep+prho[idx_xi]*pu[idx_xi]*pu[idx_xi]*pu[idx_xi])/(averStep+1.0);
                paver_rhouuv[idx_xi]=(paver_rhouuv[idx_xi]*averStep+prho[idx_xi]*pu[idx_xi]*pu[idx_xi]*pv[idx_xi])/(averStep+1.0);
                paver_rhouuw[idx_xi]=(paver_rhouuw[idx_xi]*averStep+prho[idx_xi]*pu[idx_xi]*pu[idx_xi]*pw[idx_xi])/(averStep+1.0);
                paver_rhovvu[idx_xi]=(paver_rhovvu[idx_xi]*averStep+prho[idx_xi]*pv[idx_xi]*pv[idx_xi]*pu[idx_xi])/(averStep+1.0);
                paver_rhovvv[idx_xi]=(paver_rhovvv[idx_xi]*averStep+prho[idx_xi]*pv[idx_xi]*pv[idx_xi]*pv[idx_xi])/(averStep+1.0);
                paver_rhovvw[idx_xi]=(paver_rhovvw[idx_xi]*averStep+prho[idx_xi]*pv[idx_xi]*pv[idx_xi]*pw[idx_xi])/(averStep+1.0);
                paver_rhowwu[idx_xi]=(paver_rhowwu[idx_xi]*averStep+prho[idx_xi]*pw[idx_xi]*pw[idx_xi]*pu[idx_xi])/(averStep+1.0)
                paver_rhowwv[idx_xi]=(paver_rhowwv[idx_xi]*averStep+prho[idx_xi]*pw[idx_xi]*pw[idx_xi]*pv[idx_xi])/(averStep+1.0);
                paver_rhowww[idx_xi]=(paver_rhowww[idx_xi]*averStep+prho[idx_xi]*pw[idx_xi]*pw[idx_xi]*pw[idx_xi])/(averStep+1.0);
                
                paver_dpdx[idx_xi]=(paver_dpdx[idx_xi]*averStep+ppx[idx_xi])/(averStep+1.0);
                paver_dpdy[idx_xi]=(paver_dpdy[idx_xi]*averStep+ppy[idx_xi])/(averStep+1.0);
                paver_dpdz[idx_xi]=(paver_dpdz[idx_xi]*averStep+ppz[idx_xi])/(averStep+1.0);
                
                paver_pu[idx_xi]=(paver_pu[idx_xi]*averStep+pp[idx_xi]*prhou[idx_xi])/(averStep+1.0);
                paver_pv[idx_xi]=(paver_pv[idx_xi]*averStep+pp[idx_xi]*prhov[idx_xi])/(averStep+1.0);
                paver_pw[idx_xi]=(paver_pw[idx_xi]*averStep+pp[idx_xi]*prhow[idx_xi])/(averStep+1.0);
                
                paver_dudx[idx_xi]=(paver_dudx[idx_xi]*averStep+ux[idx_xi])/(averStep+1.0);
                paver_dvdy[idx_xi]=(paver_dvdy[idx_xi]*averStep+vy[idx_xi])/(averStep+1.0);
                paver_dwdz[idx_xi]=(paver_dwdz[idx_xi]*averStep+wz[idx_xi])/(averStep+1.0);
                
                paver_pdudx[idx_xi]=(paver_pdudx[idx_xi]*averStep+pp[idx_xi]*ux[idx_xi])/(averStep+1.0);
                paver_pdvdy[idx_xi]=(paver_pdvdy[idx_xi]*averStep+pp[idx_xi]*vy[idx_xi])/(averStep+1.0);
                paver_pdwdz[idx_xi]=(paver_pdwdz[idx_xi]*averStep+pp[idx_xi]*wz[idx_xi])/(averStep+1.0);
                
                
                paver_tau_xx[idx_xi]=(paver_tau_xx[idx_xi]*averStep+ptau_xx[idx_xi])/(averStep+1.0);
                paver_tau_xy[idx_xi]=(paver_tau_xy[idx_xi]*averStep+ptau_xy[idx_xi])/(averStep+1.0);
                paver_tau_xz[idx_xi]=(paver_tau_xz[idx_xi]*averStep+ptau_xz[idx_xi])/(averStep+1.0);
                paver_tau_yy[idx_xi]=(paver_tau_yy[idx_xi]*averStep+ptau_yy[idx_xi])/(averStep+1.0);
                paver_tau_yz[idx_xi]=(paver_tau_yz[idx_xi]*averStep+ptau_yz[idx_xi])/(averStep+1.0);
                paver_tau_zz[idx_xi]=(paver_tau_zz[idx_xi]*averStep+ptau_zz[idx_xi])/(averStep+1.0);
                
            }
        }
    }
    
    averStep++;
}


void nuc3d::postproc::solveTKE(VectorField &prims,
                               VectorField &acous,
                               VectorField &xi_xyz,
                               VectorField &eta_xyz,
                               VectorField &zeta_xyz,
                               physicsModel myPhys,
                               fieldOperator3d &myOP,
                               VectorBuffer &myBf,
                               MPIComunicator3d_nonblocking &myMPI,
                               boundaryCondition &myBC)
{
    double *paver_rho=variableField[54].getDataPtr();
    double *paver_u=variableField[55].getDataPtr();
    double *paver_v=variableField[56].getDataPtr();
    double *paver_w=variableField[57].getDataPtr();
    double *paver_p=variableField[58].getDataPtr();
    double *paver_T=variableField[59].getDataPtr();
    
    //TKE budaget terms
    double *pTKE=variableField[60].getDataPtr();
    
    //
    //Production
    double *pP_tke=variableField[61].getDataPtr();
    
    double *paver_rhou=variableField[62].getDataPtr();
    double *paver_rhov=variableField[63].getDataPtr();
    double *paver_rhow=variableField[64].getDataPtr();
    
    double *paver_rhouu=variableField[65].getDataPtr();
    double *paver_rhouv=variableField[66].getDataPtr();
    double *paver_rhouw=variableField[67].getDataPtr();
    double *paver_rhovv=variableField[68].getDataPtr();
    double *paver_rhovw=variableField[69].getDataPtr();
    double *paver_rhoww=variableField[70].getDataPtr();
    
    double *pdaver_rhoudx=variableField[71].getDataPtr();
    double *pdaver_rhoudy=variableField[72].getDataPtr();
    double *pdaver_rhoudz=variableField[73].getDataPtr();
    
    double *pdaver_rhovdx=variableField[74].getDataPtr();
    double *pdaver_rhovdy=variableField[75].getDataPtr();
    double *pdaver_rhovdz=variableField[76].getDataPtr();
    
    double *pdaver_rhowdx=variableField[77].getDataPtr();
    double *pdaver_rhowdy=variableField[78].getDataPtr();
    double *pdaver_rhowdz=variableField[79].getDataPtr();
    
    //Turbulent transportation
    double *pT_tke=variableField[80].getDataPtr();
    double *pTx=variableField[81].getDataPtr();
    double *pTy=variableField[82].getDataPtr();
    double *pTz=variableField[83].getDataPtr();
    double *paver_rhouuu=variableField[84].getDataPtr();
    double *paver_rhouuv=variableField[85].getDataPtr();
    double *paver_rhouuw=variableField[86].getDataPtr();
    double *paver_rhovvu=variableField[87].getDataPtr();
    double *paver_rhovvv=variableField[88].getDataPtr();
    double *paver_rhovvw=variableField[89].getDataPtr();
    double *paver_rhowwu=variableField[90].getDataPtr();
    double *paver_rhowwv=variableField[91].getDataPtr();
    double *paver_rhowww=variableField[92].getDataPtr();
    double *pdTdx=variableField[93].getDataPtr();
    double *pdTdy=variableField[94].getDataPtr();
    double *pdTdz=variableField[95].getDataPtr();
    
    double *pPIE_tke=variableField[96].getDataPtr();
    double *pPIE_tke_0=variableField[97].getDataPtr();
    double *pPIE_tke_1=variableField[98].getDataPtr();
    double *pPIE_tke_2=variableField[99].getDataPtr();
    
    double *paver_dpdx=variableField[100].getDataPtr();
    double *paver_dpdy=variableField[101].getDataPtr();
    double *paver_dpdz=variableField[102].getDataPtr();
    
    double *paver_pu=variableField[103].getDataPtr();
    double *paver_pv=variableField[104].getDataPtr();
    double *paver_pw=variableField[105].getDataPtr();
    
    double *paver_pdudx=variableField[106].getDataPtr();
    double *paver_pdvdy=variableField[107].getDataPtr();
    double *paver_pdwdz=variableField[108].getDataPtr();
    
    double *paver_dudx=variableField[109].getDataPtr();
    double *paver_dvdy=variableField[110].getDataPtr();
    double *paver_dwdz=variableField[111].getDataPtr();
    
    double *paver_dpudx=variableField[112].getDataPtr();
    double *paver_dpvdy=variableField[113].getDataPtr();
    double *paver_dpwdz=variableField[114].getDataPtr();
    
    double *paver_pfuf=variableField[115].getDataPtr();
    double *paver_pfvf=variableField[116].getDataPtr();
    double *paver_pfwf=variableField[117].getDataPtr();
    
    double *pM_tke=variableField[118].getDataPtr();
    
    double *paver_tau_xx=variableField[119].getDataPtr();
    double *paver_tau_xy=variableField[120].getDataPtr();
    double *paver_tau_xz=variableField[121].getDataPtr();
    double *paver_tau_yy=variableField[122].getDataPtr();
    double *paver_tau_yz=variableField[123].getDataPtr();
    double *paver_tau_zz=variableField[124].getDataPtr();
    //
    double *paver_dtau_xxdx=variableField[125].getDataPtr();
    double *paver_dtau_xydy=variableField[126].getDataPtr();
    double *paver_dtau_xzdz=variableField[127].getDataPtr();
    
    double *paver_dtau_yxdx=variableField[128].getDataPtr();
    double *paver_dtau_yydy=variableField[129].getDataPtr();
    double *paver_dtau_yzdz=variableField[130].getDataPtr();
    
    double *paver_dtau_zxdx=variableField[131].getDataPtr();
    double *paver_dtau_zydy=variableField[132].getDataPtr();
    double *paver_dtau_zzdz=variableField[133].getDataPtr();
    //TKE and Production
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                
                //TKE budaget terms
                pTKE[idx_xi]=0.5*(paver_rhouu[idx_xi]/paver_rho[idx_xi]
                                  -paver_rhou[idx_xi]/paver_rho[idx_xi]*paver_rhou[idx_xi]/paver_rho[idx_xi])
                +0.5*(paver_rhovv[idx_xi]/paver_rho[idx_xi]
                      -paver_rhov[idx_xi]/paver_rho[idx_xi]*paver_rhov[idx_xi]/paver_rho[idx_xi])
                +0.5*(paver_rhoww[idx_xi]/paver_rho[idx_xi]
                      -paver_rhow[idx_xi]/paver_rho[idx_xi]*paver_rhow[idx_xi]/paver_rho[idx_xi]);
                
                //
                //Production
                pP_tke[idx_xi]=-(pdaver_rhoudx[idx_xi]*(paver_rhouu[idx_xi]/paver_rho[idx_xi]
                                                        -paver_rhou[idx_xi]/paver_rho[idx_xi]*paver_rhou[idx_xi]/paver_rho[idx_xi])
                                 +pdaver_rhoudy[idx_xi]*(paver_rhouv[idx_xi]/paver_rho[idx_xi]
                                                         -paver_rhou[idx_xi]/paver_rho[idx_xi]*paver_rhov[idx_xi]/paver_rho[idx_xi])
                                 +pdaver_rhoudz[idx_xi]*(paver_rhouw[idx_xi]/paver_rho[idx_xi]
                                                         -paver_rhou[idx_xi]/paver_rho[idx_xi]*paver_rhow[idx_xi]/paver_rho[idx_xi])
                                 +pdaver_rhovdx[idx_xi]*(paver_rhovu[idx_xi]/paver_rho[idx_xi]
                                                         -paver_rhov[idx_xi]/paver_rho[idx_xi]*paver_rhou[idx_xi]/paver_rho[idx_xi])
                                 +pdaver_rhovdy[idx_xi]*(paver_rhovv[idx_xi]/paver_rho[idx_xi]
                                                         -paver_rhov[idx_xi]/paver_rho[idx_xi]*paver_rhov[idx_xi]/paver_rho[idx_xi])
                                 +pdaver_rhovdz[idx_xi]*(paver_rhovw[idx_xi]/paver_rho[idx_xi]
                                                         -paver_rhov[idx_xi]/paver_rho[idx_xi]*paver_rhow[idx_xi]/paver_rho[idx_xi])
                                 +pdaver_rhowdx[idx_xi]*(paver_rhowu[idx_xi]/paver_rho[idx_xi]
                                                         -paver_rhow[idx_xi]/paver_rho[idx_xi]*paver_rhou[idx_xi]/paver_rho[idx_xi])
                                 +pdaver_rhowdy[idx_xi]*(paver_rhowv[idx_xi]/paver_rho[idx_xi]
                                                         -paver_rhow[idx_xi]/paver_rho[idx_xi]*paver_rhov[idx_xi]/paver_rho[idx_xi])
                                 +pdaver_rhowdz[idx_xi]*(paver_rhoww[idx_xi]/paver_rho[idx_xi]
                                                         -paver_rhow[idx_xi]/paver_rho[idx_xi]*paver_rhow[idx_xi]/paver_rho[idx_xi]));
            }
        }
    }
    
    //Turbulent Transportation
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                
                pTx[idx_xi]=0.5*(paver_rhouuu[idx_xi]+paver_rhovvu[idx_xi]+paver_rhowwu[idx_xi])
                -(paver_rhouu[idx_xi]*paver_rhou[idx_xi]/paver_rho[idx_xi]
                  +paver_rhouv[idx_xi]*paver_rhov[idx_xi]/paver_rho[idx_xi]
                  +paver_rhouw[idx_xi]*paver_rhow[idx_xi]/paver_rho[idx_xi])
                +0.5*paver_rhou[idx_xi]*(paver_rhou[idx_xi]*paver_rhou[idx_xi]/paver_rho[idx_xi]
                                         +paver_rhov[idx_xi]*paver_rhov[idx_xi]/paver_rho[idx_xi]
                                         +paver_rhow[idx_xi]*paver_rhow[idx_xi]/paver_rho[idx_xi]);
                
                pTy[idx_xi]=0.5*(paver_rhouuv[idx_xi]+paver_rhovvv[idx_xi]+paver_rhowwv[idx_xi])
                -(paver_rhouv[idx_xi]*paver_rhou[idx_xi]/paver_rho[idx_xi]
                  +paver_rhovv[idx_xi]*paver_rhov[idx_xi]/paver_rho[idx_xi]
                  +paver_rhovw[idx_xi]*paver_rhow[idx_xi]/paver_rho[idx_xi])
                +0.5*paver_rhov[idx_xi]*(paver_rhou[idx_xi]*paver_rhou[idx_xi]/paver_rho[idx_xi]
                                         +paver_rhov[idx_xi]*paver_rhov[idx_xi]/paver_rho[idx_xi]
                                         +paver_rhow[idx_xi]*paver_rhow[idx_xi]/paver_rho[idx_xi]);
                
                pTz[idx_xi]=0.5*(paver_rhouuw[idx_xi]+paver_rhovvw[idx_xi]+paver_rhowww[idx_xi])
                -(paver_rhouw[idx_xi]*paver_rhou[idx_xi]/paver_rho[idx_xi]
                  +paver_rhovw[idx_xi]*paver_rhov[idx_xi]/paver_rho[idx_xi]
                  +paver_rhoww[idx_xi]*paver_rhow[idx_xi]/paver_rho[idx_xi])
                +0.5*paver_rhow[idx_xi]*(paver_rhou[idx_xi]*paver_rhou[idx_xi]/paver_rho[idx_xi]
                                         +paver_rhov[idx_xi]*paver_rhov[idx_xi]/paver_rho[idx_xi]
                                         +paver_rhow[idx_xi]*paver_rhow[idx_xi]/paver_rho[idx_xi]);
                
                
            }
        }
    }
    
    solveGrad(variableField[81], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_x(variableField[93],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    solveGrad(variableField[82], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_y(variableField[94],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    solveGrad(variableField[83], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_z(variableField[95],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                pT_tke[idx_xi]=-(pdTdx[idx_xi]+pdTdy[idx_xi]+pdTdz[idx_xi]);
            }
        }
    }
    
    //pressure dilatation
    
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                pPIE_tke_0[idx_xi]=-(paver_u[idx_xi]-paver_rhou[idx_xi]/paver_rho[idx_xi])*paver_dpdx[idx_xi]
                -(paver_v[idx_xi]-paver_rhov[idx_xi]/paver_rho[idx_xi])*paver_dpdy[idx_xi]
                -(paver_w[idx_xi]-paver_rhow[idx_xi]/paver_rho[idx_xi])*paver_dpdz[idx_xi];
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
                paver_pfuf[idx_xi]=paver_pu[idx_xi]-paver_p[idx_xi]*paver_u[idx_xi];
                paver_pfvf[idx_xi]=paver_pv[idx_xi]-paver_p[idx_xi]*paver_v[idx_xi];
                paver_pfwf[idx_xi]=paver_pw[idx_xi]-paver_p[idx_xi]*paver_w[idx_xi];
            }
        }
    }
    
    solveGrad(variableField[115], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_x(variableField[112],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    solveGrad(variableField[116], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_y(variableField[113],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    solveGrad(variableField[117], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_z(variableField[114],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);

    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                pPIE_tke_1[idx_xi]=-paver_dpudx[idx_xi]-paver_dpvdy[idx_xi]-paver_dpwdz[idx_xi];
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
                pPIE_tke_2[idx_xi]=(paver_pdudx[idx_xi]-paver_p[idx_xi]*paver_dudx[idx_xi])
                +(paver_pdvdy[idx_xi]-paver_p[idx_xi]*paver_dvdy[idx_xi])
                +(paver_pdwdz[idx_xi]-paver_p[idx_xi]*paver_dwdz[idx_xi]);
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
                
                pPIE_tke[idx_xi]=pPIE_tke_0[idx_xi]+pPIE_tke_1[idx_xi]+pPIE_tke_2[idx_xi];
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

void nuc3d::postproc::solveGrad_x(Field Q_x,
                                  Field Q_xi,
                                  Field Q_eta,
                                  Field Q_zeta,
                                  VectorField &xi_xyz,
                                  VectorField &eta_xyz,
                                  VectorField &zeta_xyz)
{
    int nx=Q.getSizeX();
    int ny=Q.getSizeY();
    int nz=Q.getSizeZ();
    
    double *pxi_x=xi_xyz[0].getDataPtr();
    double *peta_x=eta_xyz[0].getDataPtr();
    double *pzeta_x=zeta_xyz[0].getDataPtr();
    
    double *uxi=Q_xi.getDataPtr();
    double *ueta=Q_eta.getDataPtr();
    double *uzeta=Q_zeta.getDataPtr();
    
    double *ux=Q_x.getDataPtr();
    
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
                
                double eta_x=peta_x[idx_xi];
                
                double zeta_x=pzeta_x[idx_xi];
                
                //d(u,v,w,T)/d(x,y,z)
                ux[idx_xi]=uxi[idx_xi]*xi_x+ueta[idx_eta]*eta_x+uzeta[idx_zeta]*zeta_x;
            }
        }
    }
    
    
}

void nuc3d::postproc::solveGrad_y(Field Q_y,
                                  Field Q_xi,
                                  Field Q_eta,
                                  Field Q_zeta,
                                  VectorField &xi_xyz,
                                  VectorField &eta_xyz,
                                  VectorField &zeta_xyz)
{
    
    int nx=Q.getSizeX();
    int ny=Q.getSizeY();
    int nz=Q.getSizeZ();
    
    
    double *pxi_y=xi_xyz[1].getDataPtr();
    double *peta_y=eta_xyz[1].getDataPtr();
    double *pzeta_y=zeta_xyz[1].getDataPtr();
    
    double *uxi=Q_xi.getDataPtr();
    double *ueta=Q_eta.getDataPtr();
    double *uzeta=Q_zeta.getDataPtr();
    
    double *uy=Q_y.getDataPtr();
    
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                int idx_eta=ny*nz*i+ny*k+j;
                int idx_zeta=nz*nx*j+nz*i+k;
                
                double xi_y=pxi_y[idx_xi];
                double eta_y=peta_y[idx_xi];
                double zeta_y=pzeta_y[idx_xi];
                
                //d(u,v,w,T)/d(x,y,z)
                uy[idx_xi]=uxi[idx_xi]*xi_y+ueta[idx_eta]*eta_y+uzeta[idx_zeta]*zeta_y;
                
            }
        }
    }
    
}

void nuc3d::postproc::solveGrad_z(Field Q_z,
                                  Field Q_xi,
                                  Field Q_eta,
                                  Field Q_zeta,
                                  VectorField &xi_xyz,
                                  VectorField &eta_xyz,
                                  VectorField &zeta_xyz)
{
    int nx=Q.getSizeX();
    int ny=Q.getSizeY();
    int nz=Q.getSizeZ();
    
    
    double *pxi_z=xi_xyz[2].getDataPtr();
    double *peta_z=eta_xyz[2].getDataPtr();
    double *pzeta_z=zeta_xyz[2].getDataPtr();
    
    double *uxi=Q_xi.getDataPtr();
    double *ueta=Q_eta.getDataPtr();
    double *uzeta=Q_zeta.getDataPtr();
    
    double *uz=Q_z.getDataPtr();
    
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                int idx_eta=ny*nz*i+ny*k+j;
                int idx_zeta=nz*nx*j+nz*i+k;
                
                
                double xi_z=pxi_z[idx_xi];
                double eta_z=peta_z[idx_xi];
                double zeta_z=pzeta_z[idx_xi];
                
                //d(u,v,w,T)/d(x,y,z)
                uz[idx_xi]=uxi[idx_xi]*xi_z+ueta[idx_eta]*eta_z+uzeta[idx_zeta]*zeta_z;
                
            }
        }
    }
    
}

void nuc3d::postproc::solveGrad_xyz(Field Q_x,
                                    Field Q_y,
                                    Field Q_z,
                                    Field Q_xi,
                                    Field Q_eta,
                                    Field Q_zeta,
                                    VectorField &xi_xyz,
                                    VectorField &eta_xyz,
                                    VectorField &zeta_xyz)
{
    solveGrad_x(Q_x, Q_xi, Q_eta, Q_zeta, xi_xyz, eta_xyz, zeta_xyz);
    solveGrad_y(Q_y, Q_xi, Q_eta, Q_zeta, xi_xyz, eta_xyz, zeta_xyz);
    solveGrad_z(Q_z, Q_xi, Q_eta, Q_zeta, xi_xyz, eta_xyz, zeta_xyz);
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

