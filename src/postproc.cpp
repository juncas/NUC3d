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
    /*
     PostProc.in file format
     
     0 AverageType            0:          nothing
     1:          primarys and grads
     2:          primarys, grads and derived variables
     3:          primarys, grads, derived variables and TKE
     
     1 AverageStartStep       0:          new average
     'step':     from saved file step
     
     2 SimutaneousData        0:nothing
     1:Q
     
     3 AveragedSaveType       0:binary
     1:binary and tecplot
     
     4 SimutaneousSaveType    0:not save
     1:tecplot
     
     */
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
    
    switch(variableScalarInt[0])//for average
    {
        case 0://sovle nothing
            break;
        case 1://average primarys and grads
            for(int n=0;n<VarNameList_primary.size();n++)
            {
                TemporalPrimaryField.push_back(Field(nx,ny,nz,0.0));
                AveragedPrimaryField.push_back(Field(nx,ny,nz,0.0));
            }
            
            for(int n=0;n<VarNameList_grads.size();n++)
            {
                TemporalGradsField.push_back(Field(nx,ny,nz,0.0));
                AveragedGradsField.push_back(Field(nx,ny,nz,0.0));
            }
            break;
        case 2://average primarys, grads and derived
            for(int n=0;n<VarNameList_primary.size();n++)
            {
                TemporalPrimaryField.push_back(Field(nx,ny,nz,0.0));
                AveragedPrimaryField.push_back(Field(nx,ny,nz,0.0));
            }
            
            for(int n=0;n<VarNameList_grads.size();n++)
            {
                TemporalGradsField.push_back(Field(nx,ny,nz,0.0));
                AveragedGradsField.push_back(Field(nx,ny,nz,0.0));
            }
            
            for(int n=0;n<VarNameList_derived.size();n++)
            {
                TemporalDerivedField.push_back(Field(nx,ny,nz,0.0));
                AveragedDerivedField.push_back(Field(nx,ny,nz,0.0));
            }
            
            break;
        case 3://average primarys, grads, derived and TKE
            for(int n=0;n<VarNameList_primary.size();n++)
            {
                TemporalPrimaryField.push_back(Field(nx,ny,nz,0.0));
                AveragedPrimaryField.push_back(Field(nx,ny,nz,0.0));
            }
            
            for(int n=0;n<VarNameList_grads.size();n++)
            {
                TemporalGradsField.push_back(Field(nx,ny,nz,0.0));
                AveragedGradsField.push_back(Field(nx,ny,nz,0.0));
            }
            
            for(int n=0;n<VarNameList_derived.size();n++)
            {
                TemporalDerivedField.push_back(Field(nx,ny,nz,0.0));
                AveragedDerivedField.push_back(Field(nx,ny,nz,0.0));
            }
            
            for(int n=0;n<VarNameList_TKE.size();n++)
            {
                TemporalTkeField.push_back(Field(nx,ny,nz,0.0));
                AveragedTkeField.push_back(Field(nx,ny,nz,0.0));
            }
            
            for(int n=0;n<VarNameList_TKEbudget.size();n++)
            {
                TKEbudgetField.push_back(Field(nx,ny,nz,0.0));
            }
            break;
        default:
            std::cout<<"No such postproc for "<<variableScalarInt[0]<<std::endl;
            exit(-1);
    }
    
    if(variableScalarInt[1]!=0)
    {
        std::string forename_flow = ("AveragedFlowData/Averaged_step_");
        std::string step;
        std::string mid("_ID_");
        std::string id;
        std::string tailname = (".bin");
        int myId;
        MPI_Comm_rank(MPI_COMM_WORLD, &myId);
        
        std::stringstream ss_step,ss_id;
        ss_step << variableScalarInt[3];
        ss_step >> step;
        
        ss_id<<myId;
        ss_id>>id;
        
        std::string filename_flow = forename_flow + step + mid + id + tailname;
        
        std::ifstream myIOfile;
        
        myIOfile.open(filename_flow,std::ios::in|std::ios::binary);
        myIOfile.read(reinterpret_cast<char *>(&averStep), sizeof (averStep));
        
        switch(variableScalarInt[0])//for average
        {
            case 0:
                break;
            case 1://averaged primarys and grads
                for(int n=0;n<VarNameList_primary.size();n++)
                {
                    Field &data=AveragedPrimaryField[n];
                    readField_binary(myIOfile, data);
                }
                
                for(int n=0;n<VarNameList_grads.size();n++)
                {
                    Field &data=AveragedGradsField[n];
                    readField_binary(myIOfile, data);
                }
                break;
            case 2://averaged primarys, grads and derived
                for(int n=0;n<VarNameList_primary.size();n++)
                {
                    Field &data=AveragedPrimaryField[n];
                    readField_binary(myIOfile, data);
                }
                
                for(int n=0;n<VarNameList_grads.size();n++)
                {
                    Field &data=AveragedGradsField[n];
                    readField_binary(myIOfile, data);
                }
                
                for(int n=0;n<VarNameList_derived.size();n++)
                {
                    Field &data=AveragedDerivedField[n];
                    readField_binary(myIOfile, data);
                }
                break;
            case 3://average primarys, grads, derived and TKE
                for(int n=0;n<VarNameList_primary.size();n++)
                {
                    Field &data=AveragedPrimaryField[n];
                    readField_binary(myIOfile, data);
                }
                
                for(int n=0;n<VarNameList_grads.size();n++)
                {
                    Field &data=AveragedGradsField[n];
                    readField_binary(myIOfile, data);
                }
                
                for(int n=0;n<VarNameList_derived.size();n++)
                {
                    Field &data=AveragedDerivedField[n];
                    readField_binary(myIOfile, data);
                }
                
                for(int n=0;n<VarNameList_TKE.size();n++)
                {
                    Field &data=AveragedTkeField[n];
                    readField_binary(myIOfile, data);
                }
                
                for(int n=0;n<VarNameList_TKEbudget.size();n++)
                {
                    Field &data=TKEbudgetField[n];
                    readField_binary(myIOfile, data);
                }
                break;
            default:
                std::cout<<"No such postproc data for "<<variableScalarInt[0]<<std::endl;
                exit(-1);
        }
    }
    else
    {
        averStep=0;
    }
    
    switch(variableScalarInt[2])//for simutaneous variables
    {
        case 0:// do nonthing
            break;
        case 1://sovle Q
            for(int n=0;n<VarNameList_q.size();n++)
            {
                TemporalQField.push_back(Field(nx,ny,nz,0.0));
            }
            break;
        default:
            std::cout<<"No such postproc for "<<variableScalarInt[1]<<std::endl;
            exit(-1);
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
                                physicsModel &myPhys,
                                fieldOperator3d &myOP,
                                VectorBuffer &myBf,
                                MPIComunicator3d_nonblocking &myMPI,
                                boundaryCondition &myBC,
                                IOController &myIO,
                                int istep,
                                double time)
{
    if(variableScalarInt[0]!=0)
    {
        solveAveraged(prims, acous,xi_xyz, eta_xyz, zeta_xyz,myPhys, myOP, myBf, myMPI, myBC);
    }
    
    if(variableScalarInt[2]!=0)
    {
        solveTemporal(prims, acous,xi_xyz, eta_xyz, zeta_xyz,myPhys, myOP, myBf, myMPI, myBC,istep,time);
    }
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
    if(variableScalarInt[0]!=0)
    {
        switch(variableScalarInt[3])//for averaged variables
        {
            case 0:
                OutputAveraged_bin(prims, acous, xyz, xi_xyz, eta_xyz, zeta_xyz, myOP, myBf, myMPI, myBC, myIO, istep, time);
                break;
            case 1:
                OutputAveraged_bin(prims, acous, xyz, xi_xyz, eta_xyz, zeta_xyz, myOP, myBf, myMPI, myBC, myIO, istep, time);
                OutputAveraged_tec(prims, acous, xyz, xi_xyz, eta_xyz, zeta_xyz, myOP, myBf, myMPI, myBC, myIO, istep, time);
                break;
            default:
                OutputAveraged_bin(prims, acous, xyz, xi_xyz, eta_xyz, zeta_xyz, myOP, myBf, myMPI, myBC, myIO, istep, time);
                break;
        }
    }
    if(variableScalarInt[4]!=0)//for simutaneous variables
    {
        OutputTemporal(prims, acous, xyz, xi_xyz, eta_xyz, zeta_xyz, myOP, myBf, myMPI, myBC, myIO, istep, time);
    }
    
}


void nuc3d::postproc::OutputTemporal(VectorField &prims,
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
    std::string forename_flow = ("flowDataTemp/Q_tec_step_");
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
    
    for(int i=0;i<VarNameList_q.size();i++)
    {
        myIOfile<<","<<VarNameList_q[i];
    }
    
    myIOfile<<"\n Zone I = "<<nx+1<<", J= "<<ny+1<<", K="<<nz+1
    <<"\n DATAPACKING=BLOCK, VARLOCATION=(["<<xyz.size()+1<<"-"
    <<xyz.size()+VarNameList_q.size()
    <<"]=CELLCENTERED)\n";
    
    for(auto iter=xyz.begin();iter!=xyz.end();iter++)
    {
        writeField(myIOfile, *iter);
    }
    
    for(int n=0;n<VarNameList_q.size();n++)
    {
        writeField(myIOfile,TemporalQField[n]);
    }
    
    /* ADD arbitrary number of post values
     writeField(myIOfile, omega2);
     */
    
    
    myIOfile.close();
    
    
}

void nuc3d::postproc::OutputAveraged_bin(VectorField &prims,
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
    std::string forename_flow = ("AveragedFlowData/Averaged_step_");
    std::string step;
    std::string mid("_id_");
    std::string id;
    std::string tailname = (".bin");
    int myID=myMPI.getMyId();
    
    std::stringstream ss_step,ss_id;
    ss_step << istep;
    ss_step >> step;
    
    ss_id<<myID;
    ss_id>>id;
    
    std::string filename_flow = forename_flow + step + mid + id + tailname;
    
    std::ofstream myIOfile;
    
    myIOfile.open(filename_flow,std::ios::out|std::ios::binary);
    
    myIOfile.write(reinterpret_cast<char *>(&averStep), sizeof(averStep));
    
    switch(variableScalarInt[0])//for average
    {
            
        case 1://averaged primarys and grads
            for(int n=0;n<VarNameList_primary.size();n++)
            {
                Field &data=AveragedPrimaryField[n];
                writeField_binary(myIOfile, data);
            }
            
            for(int n=0;n<VarNameList_grads.size();n++)
            {
                Field &data=AveragedGradsField[n];
                writeField_binary(myIOfile, data);
            }
            break;
        case 2://averaged primarys, grads and derived
            for(int n=0;n<VarNameList_primary.size();n++)
            {
                Field &data=AveragedPrimaryField[n];
                writeField_binary(myIOfile, data);
            }
            
            for(int n=0;n<VarNameList_grads.size();n++)
            {
                Field &data=AveragedGradsField[n];
                writeField_binary(myIOfile, data);
            }
            
            for(int n=0;n<VarNameList_derived.size();n++)
            {
                Field &data=AveragedDerivedField[n];
                writeField_binary(myIOfile, data);
            }
            break;
        case 3://average primarys, grads, derived and TKE
            for(int n=0;n<VarNameList_primary.size();n++)
            {
                Field &data=AveragedPrimaryField[n];
                writeField_binary(myIOfile, data);
            }
            
            for(int n=0;n<VarNameList_grads.size();n++)
            {
                Field &data=AveragedGradsField[n];
                writeField_binary(myIOfile, data);
            }
            
            for(int n=0;n<VarNameList_derived.size();n++)
            {
                Field &data=AveragedDerivedField[n];
                writeField_binary(myIOfile, data);
            }
            
            for(int n=0;n<VarNameList_TKE.size();n++)
            {
                Field &data=AveragedTkeField[n];
                writeField_binary(myIOfile, data);
            }
            
            for(int n=0;n<VarNameList_TKEbudget.size();n++)
            {
                Field &data=TKEbudgetField[n];
                writeField_binary(myIOfile, data);
            }
            break;
        default:
            std::cout<<"No such postproc data for "<<variableScalarInt[0]<<std::endl;
            exit(-1);
    }
    
    myIOfile.close();
    
    /* ADD arbitrary number of post values
     writeField(myIOfile, omega2);
     */
    
}

void nuc3d::postproc::OutputAveraged_tec(VectorField &prims,
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
    
    std::string forename_flow = ("AveragedFlowData/Averaged_tec_step_");
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
    
    unsigned long varnum = 0;
    
    switch(variableScalarInt[0])//for average
    {
        case 1:
            varnum=xyz.size()
            +VarNameList_primary.size()
            +VarNameList_grads.size();
            
            for(int i=0;i<VarNameList_primary.size();i++)
            {
                myIOfile<<","<<VarNameList_primary[i];
            }
            
            for(int i=0;i<VarNameList_grads.size();i++)
            {
                myIOfile<<","<<VarNameList_grads[i];
            }

            break;
        case 2:
            varnum=xyz.size()
            +VarNameList_primary.size()
            +VarNameList_grads.size()
            +VarNameList_derived.size();
            
            
            for(int i=0;i<VarNameList_primary.size();i++)
            {
                myIOfile<<","<<VarNameList_primary[i];
            }
            
            for(int i=0;i<VarNameList_grads.size();i++)
            {
                myIOfile<<","<<VarNameList_grads[i];
            }
            
            for(int i=0;i<VarNameList_derived.size();i++)
            {
                myIOfile<<","<<VarNameList_derived[i];
            }

            break;
        case 3:
            varnum=xyz.size()
            +VarNameList_primary.size()
            +VarNameList_grads.size()
            +VarNameList_derived.size()
            +VarNameList_TKE.size()
            +VarNameList_TKEbudget.size();
            
            for(int i=0;i<VarNameList_primary.size();i++)
            {
                myIOfile<<","<<VarNameList_primary[i];
            }
            
            for(int i=0;i<VarNameList_grads.size();i++)
            {
                myIOfile<<","<<VarNameList_grads[i];
            }
            
            for(int i=0;i<VarNameList_derived.size();i++)
            {
                myIOfile<<","<<VarNameList_derived[i];
            }
            
            for(int i=0;i<VarNameList_TKE.size();i++)
            {
                myIOfile<<","<<VarNameList_TKE[i];
            }
            
            for(int i=0;i<VarNameList_TKEbudget.size();i++)
            {
                myIOfile<<","<<VarNameList_TKEbudget[i];
            }

            break;
    }
    myIOfile<<"\n Zone I = "<<nx+1<<", J= "<<ny+1<<", K="<<nz+1
    <<"\n DATAPACKING=BLOCK, VARLOCATION=(["<<xyz.size()+1<<"-"
    <<varnum
    <<"]=CELLCENTERED)\n";
    
    for(auto iter=xyz.begin();iter!=xyz.end();iter++)
    {
        writeField(myIOfile, *iter);
    }
    
    switch(variableScalarInt[0])//for average
    {
        case 1://averaged primarys and grads
            for(int n=0;n<VarNameList_primary.size();n++)
            {
                Field &data=AveragedPrimaryField[n];
                writeField(myIOfile, data);
            }
            
            for(int n=0;n<VarNameList_grads.size();n++)
            {
                Field &data=AveragedGradsField[n];
                writeField(myIOfile, data);
            }
            break;
        case 2://averaged primarys, grads and derived
            for(int n=0;n<VarNameList_primary.size();n++)
            {
                Field &data=AveragedPrimaryField[n];
                writeField(myIOfile, data);
            }
            
            for(int n=0;n<VarNameList_grads.size();n++)
            {
                Field &data=AveragedGradsField[n];
                writeField(myIOfile, data);
            }
            
            for(int n=0;n<VarNameList_derived.size();n++)
            {
                Field &data=AveragedDerivedField[n];
                writeField(myIOfile, data);
            }
            break;
        case 3://average primarys, grads, derived and TKE
            for(int n=0;n<VarNameList_primary.size();n++)
            {
                Field &data=AveragedPrimaryField[n];
                writeField(myIOfile, data);
            }
            
            for(int n=0;n<VarNameList_grads.size();n++)
            {
                Field &data=AveragedGradsField[n];
                writeField(myIOfile, data);
            }
            
            for(int n=0;n<VarNameList_derived.size();n++)
            {
                Field &data=AveragedDerivedField[n];
                writeField(myIOfile, data);
            }
            
            for(int n=0;n<VarNameList_TKE.size();n++)
            {
                Field &data=AveragedTkeField[n];
                writeField(myIOfile, data);
            }
            
            for(int n=0;n<VarNameList_TKEbudget.size();n++)
            {
                Field &data=TKEbudgetField[n];
                writeField(myIOfile, data);
            }
            break;
        default:
            std::cout<<"No such postproc data for "<<variableScalarInt[0]<<std::endl;
            exit(-1);
    }
    
    myIOfile.close();
}

void nuc3d::postproc::solveTemporal(VectorField &prims,
                                    VectorField &acous,
                                    VectorField &xi_xyz,
                                    VectorField &eta_xyz,
                                    VectorField &zeta_xyz,
                                    physicsModel &myPhys,
                                    fieldOperator3d &myOP,
                                    VectorBuffer &myBf,
                                    MPIComunicator3d_nonblocking &myMPI,
                                    boundaryCondition &myBC,
                                    int istep,
                                    double time)
{
    solveQ(prims,xi_xyz,eta_xyz,zeta_xyz,myOP,myBf,myMPI,myBC,istep,time);
    
}

void nuc3d::postproc::solveQ(VectorField &prims,
                             VectorField &xi_xyz,
                             VectorField &eta_xyz,
                             VectorField &zeta_xyz,
                             fieldOperator3d &myOP,
                             VectorBuffer &myBf,
                             MPIComunicator3d_nonblocking &myMPI,
                             boundaryCondition &myBC,
                             int istep,
                             double time)
{
    Field &rho=prims[0];
    Field &u=prims[1];
    Field &v=prims[2];
    Field &w=prims[3];
    
    solveGrad(u, myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(TemporalQField[0],TemporalQField[1],TemporalQField[2],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    solveGrad(v, myOP, myBf[1], myMPI,myBC,1,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(TemporalQField[3],TemporalQField[4],TemporalQField[5],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    solveGrad(w, myOP, myBf[2], myMPI,myBC,2,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(TemporalQField[6],TemporalQField[7],TemporalQField[8],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    double *prho=rho.getDataPtr();
    double *pu=u.getDataPtr();
    double *pv=v.getDataPtr();
    double *pw=w.getDataPtr();
    
    double *ux=TemporalQField[0].getDataPtr();
    double *uy=TemporalQField[1].getDataPtr();
    double *uz=TemporalQField[2].getDataPtr();
    
    double *vx=TemporalQField[3].getDataPtr();
    double *vy=TemporalQField[4].getDataPtr();
    double *vz=TemporalQField[5].getDataPtr();
    
    double *wx=TemporalQField[6].getDataPtr();
    double *wy=TemporalQField[7].getDataPtr();
    double *wz=TemporalQField[8].getDataPtr();
    
    double *pOmega0=TemporalQField[9].getDataPtr();
    double *pOmega1=TemporalQField[10].getDataPtr();
    double *pOmega2=TemporalQField[11].getDataPtr();
    double *pEns=TemporalQField[12].getDataPtr();
    double *pQ=TemporalQField[13].getDataPtr();
    
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
                pEns[idx_xi]=w0*w0+w1*w1+w2*w2;
                pQ[idx_xi]=q0;
                enstrophy+=0.5*(w0*w0+w1*w1+w2*w2);
                kinetic+=0.5*prho[idx_xi]*(pu[idx_xi]*pu[idx_xi]
                                           +pv[idx_xi]*pv[idx_xi]
                                           +pw[idx_xi]*pw[idx_xi]);
            }
        }
    }
    
    MPI_Allreduce(&enstrophy, &enstrophy_glb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    MPI_Allreduce(&kinetic, &kinetic_glb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    if(0==myMPI.getMyId())
    {
        std::ofstream myIOfile;
        myIOfile.open("enstrophy.dat",std::ios::out|std::ios::app);
        myIOfile<<std::setprecision(12)
        <<istep
        <<" "
        <<time
        <<" "
        <<enstrophy_glb
        <<" "
        <<kinetic_glb<<"\n";
        myIOfile.close();
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
}
//---------------------------------------------------
void nuc3d::postproc::solveAveraged(VectorField &prims,
                                    VectorField &acous,
                                    VectorField &xi_xyz,
                                    VectorField &eta_xyz,
                                    VectorField &zeta_xyz,
                                    physicsModel &myPhys,
                                    fieldOperator3d &myOP,
                                    VectorBuffer &myBf,
                                    MPIComunicator3d_nonblocking &myMPI,
                                    boundaryCondition &myBC)
{
    for(int n=0;n<variableScalarInt[0];n++)
    {
        (this->*myAverage[n])(prims,
                              acous,
                              xi_xyz,
                              eta_xyz,
                              zeta_xyz,
                              myPhys,
                              myOP,
                              myBf,
                              myMPI,
                              myBC);
    }
    
    averStep++;
}

void nuc3d::postproc::solveAveraged0(VectorField &prims,
                                     VectorField &acous,
                                     VectorField &xi_xyz,
                                     VectorField &eta_xyz,
                                     VectorField &zeta_xyz,
                                     physicsModel &myPhys,
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
    Field &e=acous[1];
    Field &alpha=acous[2];
    
    myPhys.getMiu(T, TemporalPrimaryField[8], TemporalPrimaryField[9]);
    
    solveGrad(u, myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(TemporalGradsField[0],TemporalGradsField[1],TemporalGradsField[2],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    solveGrad(v, myOP, myBf[1], myMPI,myBC,1,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(TemporalGradsField[3],TemporalGradsField[4],TemporalGradsField[5],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    solveGrad(w, myOP, myBf[2], myMPI,myBC,2,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(TemporalGradsField[6],TemporalGradsField[7],TemporalGradsField[8],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    solveGrad(p, myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(TemporalGradsField[9],TemporalGradsField[10],TemporalGradsField[11],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    solveGrad(T, myOP, myBf[1], myMPI,myBC,1,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(TemporalGradsField[12],TemporalGradsField[13],TemporalGradsField[14],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    
    double *prho=rho.getDataPtr();
    double *pu=u.getDataPtr();
    double *pv=v.getDataPtr();
    double *pw=w.getDataPtr();
    double *pp=p.getDataPtr();
    double *pT=T.getDataPtr();
    double *pe=e.getDataPtr();
    double *palpha=alpha.getDataPtr();
    
    double* ux_temp=TemporalGradsField[0].getDataPtr();
    double* uy_temp=TemporalGradsField[1].getDataPtr();
    double* uz_temp=TemporalGradsField[2].getDataPtr();
    double* vx_temp=TemporalGradsField[3].getDataPtr();
    double* vy_temp=TemporalGradsField[4].getDataPtr();
    double* vz_temp=TemporalGradsField[5].getDataPtr();
    double* wx_temp=TemporalGradsField[6].getDataPtr();
    double* wy_temp=TemporalGradsField[7].getDataPtr();
    double* wz_temp=TemporalGradsField[8].getDataPtr();
    double *px_temp=TemporalGradsField[9].getDataPtr();
    double *py_temp=TemporalGradsField[10].getDataPtr();
    double *pz_temp=TemporalGradsField[11].getDataPtr();
    double *tx_temp=TemporalGradsField[12].getDataPtr();
    double *ty_temp=TemporalGradsField[13].getDataPtr();
    double *tz_temp=TemporalGradsField[14].getDataPtr();
    double *Omega0_temp=TemporalGradsField[15].getDataPtr();
    double *Omega1_temp=TemporalGradsField[16].getDataPtr();
    double *Omega2_temp=TemporalGradsField[17].getDataPtr();
    
    double *rho_temp=TemporalPrimaryField[0].getDataPtr();
    double *u_temp=TemporalPrimaryField[1].getDataPtr();
    double *v_temp=TemporalPrimaryField[2].getDataPtr();
    double *w_temp=TemporalPrimaryField[3].getDataPtr();
    double *p_temp=TemporalPrimaryField[4].getDataPtr();
    double *T_temp=TemporalPrimaryField[5].getDataPtr();
    double *e_temp=TemporalPrimaryField[6].getDataPtr();
    double *alpha_temp=TemporalPrimaryField[7].getDataPtr();
    double* miu_temp=TemporalPrimaryField[8].getDataPtr();
    
    
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                
                rho_temp[idx_xi]=prho[idx_xi];
                u_temp[idx_xi]=pu[idx_xi];
                v_temp[idx_xi]=pv[idx_xi];
                w_temp[idx_xi]=pw[idx_xi];
                p_temp[idx_xi]=pp[idx_xi];
                T_temp[idx_xi]=pT[idx_xi];
                e_temp[idx_xi]=pe[idx_xi];
                alpha_temp[idx_xi]=palpha[idx_xi];
                
                Omega0_temp[idx_xi]=wy_temp[idx_xi]-vz_temp[idx_xi];
                Omega1_temp[idx_xi]=uz_temp[idx_xi]-wx_temp[idx_xi];
                Omega2_temp[idx_xi]=vx_temp[idx_xi]-uy_temp[idx_xi];
            }
        }
    }
    
    
    for(int n=0;n<VarNameList_primary.size();n++)
    {
        double *f_aver=AveragedPrimaryField[n].getDataPtr();
        double *f_temp=TemporalPrimaryField[n].getDataPtr();
        
        for (int k=0; k<nz; k++)
        {
            for (int j=0; j<ny; j++)
            {
                for (int i=0; i<nx; i++)
                {
                    int idx_xi=nx*ny*k+nx*j+i;
                    
                    f_aver[idx_xi]=(f_aver[idx_xi]*averStep+f_temp[idx_xi])/(averStep+1.0);
                }
            }
        }
    }
    
    for(int n=0;n<VarNameList_grads.size();n++)
    {
        double *f_aver=AveragedGradsField[n].getDataPtr();
        double *f_temp=TemporalGradsField[n].getDataPtr();
        
        for (int k=0; k<nz; k++)
        {
            for (int j=0; j<ny; j++)
            {
                for (int i=0; i<nx; i++)
                {
                    int idx_xi=nx*ny*k+nx*j+i;
                    
                    f_aver[idx_xi]=(f_aver[idx_xi]*averStep+f_temp[idx_xi])/(averStep+1.0);
                }
            }
        }
    }
}

void nuc3d::postproc::solveAveraged1(VectorField &prims,
                                     VectorField &acous,
                                     VectorField &xi_xyz,
                                     VectorField &eta_xyz,
                                     VectorField &zeta_xyz,
                                     physicsModel &myPhys,
                                     fieldOperator3d &myOP,
                                     VectorBuffer &myBf,
                                     MPIComunicator3d_nonblocking &myMPI,
                                     boundaryCondition &myBC)
{
    double *rho_temp=TemporalPrimaryField[0].getDataPtr();
    double *u_temp=TemporalPrimaryField[1].getDataPtr();
    double *v_temp=TemporalPrimaryField[2].getDataPtr();
    double *w_temp=TemporalPrimaryField[3].getDataPtr();
    double *p_temp=TemporalPrimaryField[4].getDataPtr();
    double *T_temp=TemporalPrimaryField[5].getDataPtr();
    double *e_temp=TemporalPrimaryField[6].getDataPtr();
    double *alpha_temp=TemporalPrimaryField[7].getDataPtr();
    double* miu_temp=TemporalPrimaryField[8].getDataPtr();
    
    double* ux_temp=TemporalGradsField[0].getDataPtr();
    double* uy_temp=TemporalGradsField[1].getDataPtr();
    double* uz_temp=TemporalGradsField[2].getDataPtr();
    double* vx_temp=TemporalGradsField[3].getDataPtr();
    double* vy_temp=TemporalGradsField[4].getDataPtr();
    double* vz_temp=TemporalGradsField[5].getDataPtr();
    double* wx_temp=TemporalGradsField[6].getDataPtr();
    double* wy_temp=TemporalGradsField[7].getDataPtr();
    double* wz_temp=TemporalGradsField[8].getDataPtr();
    double *px_temp=TemporalGradsField[9].getDataPtr();
    double *py_temp=TemporalGradsField[10].getDataPtr();
    double *pz_temp=TemporalGradsField[11].getDataPtr();
    double *tx_temp=TemporalGradsField[12].getDataPtr();
    double *ty_temp=TemporalGradsField[13].getDataPtr();
    double *tz_temp=TemporalGradsField[14].getDataPtr();
    double *Omega0_temp=TemporalGradsField[15].getDataPtr();
    double *Omega1_temp=TemporalGradsField[16].getDataPtr();
    double *Omega2_temp=TemporalGradsField[17].getDataPtr();
    
    double *rhou_temp=TemporalDerivedField[0].getDataPtr();
    double *rhov_temp=TemporalDerivedField[1].getDataPtr();
    double *rhow_temp=TemporalDerivedField[2].getDataPtr();
    double *rhoe_temp=TemporalDerivedField[3].getDataPtr();
    double *uu_temp=TemporalDerivedField[4].getDataPtr();
    double *uv_temp=TemporalDerivedField[5].getDataPtr();
    double *uw_temp=TemporalDerivedField[6].getDataPtr();
    double *vv_temp=TemporalDerivedField[7].getDataPtr();
    double *vw_temp=TemporalDerivedField[8].getDataPtr();
    double *ww_temp=TemporalDerivedField[9].getDataPtr();
    double *rhouu_temp=TemporalDerivedField[10].getDataPtr();
    double *rhouv_temp=TemporalDerivedField[11].getDataPtr();
    double *rhouw_temp=TemporalDerivedField[12].getDataPtr();
    double *rhovv_temp=TemporalDerivedField[13].getDataPtr();
    double *rhovw_temp=TemporalDerivedField[14].getDataPtr();
    double *rhoww_temp=TemporalDerivedField[15].getDataPtr();
    double *tau_xx_temp=TemporalDerivedField[16].getDataPtr();
    double *tau_xy_temp=TemporalDerivedField[17].getDataPtr();
    double *tau_xz_temp=TemporalDerivedField[18].getDataPtr();
    double *tau_yy_temp=TemporalDerivedField[19].getDataPtr();
    double *tau_yz_temp=TemporalDerivedField[20].getDataPtr();
    double *tau_zz_temp=TemporalDerivedField[21].getDataPtr();
    double *rhotau_xx_temp=TemporalDerivedField[22].getDataPtr();
    double *rhotau_xy_temp=TemporalDerivedField[23].getDataPtr();
    double *rhotau_xz_temp=TemporalDerivedField[24].getDataPtr();
    double *rhotau_yy_temp=TemporalDerivedField[25].getDataPtr();
    double *rhotau_yz_temp=TemporalDerivedField[26].getDataPtr();
    double *rhotau_zz_temp=TemporalDerivedField[27].getDataPtr();
    
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                
                rhou_temp[idx_xi]=rho_temp[idx_xi]*u_temp[idx_xi];
                rhov_temp[idx_xi]=rho_temp[idx_xi]*v_temp[idx_xi];
                rhow_temp[idx_xi]=rho_temp[idx_xi]*w_temp[idx_xi];
                rhoe_temp[idx_xi]=rho_temp[idx_xi]*e_temp[idx_xi];
                uu_temp[idx_xi]=u_temp[idx_xi]*u_temp[idx_xi];
                uv_temp[idx_xi]=u_temp[idx_xi]*v_temp[idx_xi];
                uw_temp[idx_xi]=u_temp[idx_xi]*w_temp[idx_xi];
                vv_temp[idx_xi]=v_temp[idx_xi]*v_temp[idx_xi];
                vw_temp[idx_xi]=v_temp[idx_xi]*w_temp[idx_xi];
                ww_temp[idx_xi]=w_temp[idx_xi]*w_temp[idx_xi];
                
                rhouu_temp[idx_xi]=rho_temp[idx_xi]*u_temp[idx_xi]*u_temp[idx_xi];
                rhouv_temp[idx_xi]=rho_temp[idx_xi]*u_temp[idx_xi]*u_temp[idx_xi];
                rhouw_temp[idx_xi]=rho_temp[idx_xi]*u_temp[idx_xi]*u_temp[idx_xi];
                rhovv_temp[idx_xi]=rho_temp[idx_xi]*v_temp[idx_xi]*v_temp[idx_xi];
                rhovw_temp[idx_xi]=rho_temp[idx_xi]*v_temp[idx_xi]*v_temp[idx_xi];
                rhoww_temp[idx_xi]=rho_temp[idx_xi]*w_temp[idx_xi]*w_temp[idx_xi];
                
                double velgrad=ux_temp[idx_xi]+vy_temp[idx_xi]+wz_temp[idx_xi];
                
                tau_xx_temp[idx_xi]=miu_temp[idx_xi]*(2.0*ux_temp[idx_xi]-2.0/3.0*velgrad);
                tau_xy_temp[idx_xi]=miu_temp[idx_xi]*(uy_temp[idx_xi]+vx_temp[idx_xi]-2.0/3.0*velgrad);
                tau_xz_temp[idx_xi]=miu_temp[idx_xi]*(uz_temp[idx_xi]+wx_temp[idx_xi]-2.0/3.0*velgrad);
                tau_yy_temp[idx_xi]=miu_temp[idx_xi]*(2.0*vy_temp[idx_xi]-2.0/3.0*velgrad);
                tau_yz_temp[idx_xi]=miu_temp[idx_xi]*(vz_temp[idx_xi]+wy_temp[idx_xi]-2.0/3.0*velgrad);
                tau_zz_temp[idx_xi]=miu_temp[idx_xi]*(2.0*wz_temp[idx_xi]-2.0/3.0*velgrad);
                
                rhotau_xx_temp[idx_xi]=rho_temp[idx_xi]*tau_xx_temp[idx_xi];
                rhotau_xy_temp[idx_xi]=rho_temp[idx_xi]*tau_xy_temp[idx_xi];
                rhotau_xz_temp[idx_xi]=rho_temp[idx_xi]*tau_xz_temp[idx_xi];
                rhotau_yy_temp[idx_xi]=rho_temp[idx_xi]*tau_yy_temp[idx_xi];
                rhotau_yz_temp[idx_xi]=rho_temp[idx_xi]*tau_yz_temp[idx_xi];
                rhotau_zz_temp[idx_xi]=rho_temp[idx_xi]*tau_zz_temp[idx_xi];
                
            }
        }
    }
    
    for(int n=0;n<VarNameList_derived.size();n++)
    {
        double *f_aver=AveragedDerivedField[n].getDataPtr();
        double *f_temp=TemporalDerivedField[n].getDataPtr();
        
        for (int k=0; k<nz; k++)
        {
            for (int j=0; j<ny; j++)
            {
                for (int i=0; i<nx; i++)
                {
                    int idx_xi=nx*ny*k+nx*j+i;
                    
                    f_aver[idx_xi]=(f_aver[idx_xi]*averStep+f_temp[idx_xi])/(averStep+1.0);
                }
            }
        }
    }
}

void nuc3d::postproc::solveAveraged2(VectorField &prims,
                                     VectorField &acous,
                                     VectorField &xi_xyz,
                                     VectorField &eta_xyz,
                                     VectorField &zeta_xyz,
                                     physicsModel &myPhys,
                                     fieldOperator3d &myOP,
                                     VectorBuffer &myBf,
                                     MPIComunicator3d_nonblocking &myMPI,
                                     boundaryCondition &myBC)
{
    double *rho_temp=TemporalPrimaryField[0].getDataPtr();
    double *u_temp=TemporalPrimaryField[1].getDataPtr();
    double *v_temp=TemporalPrimaryField[2].getDataPtr();
    double *w_temp=TemporalPrimaryField[3].getDataPtr();
    double *p_temp=TemporalPrimaryField[4].getDataPtr();
    double *T_temp=TemporalPrimaryField[5].getDataPtr();
    double *e_temp=TemporalPrimaryField[6].getDataPtr();
    double *alpha_temp=TemporalPrimaryField[7].getDataPtr();
    double* miu_temp=TemporalPrimaryField[8].getDataPtr();
    
    double* ux_temp=TemporalGradsField[0].getDataPtr();
    double* uy_temp=TemporalGradsField[1].getDataPtr();
    double* uz_temp=TemporalGradsField[2].getDataPtr();
    double* vx_temp=TemporalGradsField[3].getDataPtr();
    double* vy_temp=TemporalGradsField[4].getDataPtr();
    double* vz_temp=TemporalGradsField[5].getDataPtr();
    double* wx_temp=TemporalGradsField[6].getDataPtr();
    double* wy_temp=TemporalGradsField[7].getDataPtr();
    double* wz_temp=TemporalGradsField[8].getDataPtr();
    double *px_temp=TemporalGradsField[9].getDataPtr();
    double *py_temp=TemporalGradsField[10].getDataPtr();
    double *pz_temp=TemporalGradsField[11].getDataPtr();
    double *tx_temp=TemporalGradsField[12].getDataPtr();
    double *ty_temp=TemporalGradsField[13].getDataPtr();
    double *tz_temp=TemporalGradsField[14].getDataPtr();
    double *Omega0_temp=TemporalGradsField[15].getDataPtr();
    double *Omega1_temp=TemporalGradsField[16].getDataPtr();
    double *Omega2_temp=TemporalGradsField[17].getDataPtr();
    
    double *rhou_temp=TemporalDerivedField[0].getDataPtr();
    double *rhov_temp=TemporalDerivedField[1].getDataPtr();
    double *rhow_temp=TemporalDerivedField[2].getDataPtr();
    double *rhoe_temp=TemporalDerivedField[3].getDataPtr();
    double *uu_temp=TemporalDerivedField[4].getDataPtr();
    double *uv_temp=TemporalDerivedField[5].getDataPtr();
    double *uw_temp=TemporalDerivedField[6].getDataPtr();
    double *vv_temp=TemporalDerivedField[7].getDataPtr();
    double *vw_temp=TemporalDerivedField[8].getDataPtr();
    double *ww_temp=TemporalDerivedField[9].getDataPtr();
    double *rhouu_temp=TemporalDerivedField[10].getDataPtr();
    double *rhouv_temp=TemporalDerivedField[11].getDataPtr();
    double *rhouw_temp=TemporalDerivedField[12].getDataPtr();
    double *rhovv_temp=TemporalDerivedField[13].getDataPtr();
    double *rhovw_temp=TemporalDerivedField[14].getDataPtr();
    double *rhoww_temp=TemporalDerivedField[15].getDataPtr();
    double *tau_xx_temp=TemporalDerivedField[16].getDataPtr();
    double *tau_xy_temp=TemporalDerivedField[17].getDataPtr();
    double *tau_xz_temp=TemporalDerivedField[18].getDataPtr();
    double *tau_yy_temp=TemporalDerivedField[19].getDataPtr();
    double *tau_yz_temp=TemporalDerivedField[20].getDataPtr();
    double *tau_zz_temp=TemporalDerivedField[21].getDataPtr();
    double *rhotau_xx_temp=TemporalDerivedField[22].getDataPtr();
    double *rhotau_xy_temp=TemporalDerivedField[23].getDataPtr();
    double *rhotau_xz_temp=TemporalDerivedField[24].getDataPtr();
    double *rhotau_yy_temp=TemporalDerivedField[25].getDataPtr();
    double *rhotau_yz_temp=TemporalDerivedField[26].getDataPtr();
    double *rhotau_zz_temp=TemporalDerivedField[27].getDataPtr();
    
    double *up_temp=TemporalTkeField[0].getDataPtr();
    double *vp_temp=TemporalTkeField[1].getDataPtr();
    double *wp_temp=TemporalTkeField[2].getDataPtr();
    double *pdudx_temp=TemporalTkeField[3].getDataPtr();
    double *pdvdy_temp=TemporalTkeField[4].getDataPtr();
    double *pdwdz_temp=TemporalTkeField[5].getDataPtr();
    
    double *dtau_xxdx_temp=TemporalTkeField[6].getDataPtr();
    double *dtau_yxdx_temp=TemporalTkeField[7].getDataPtr();
    double *dtau_zxdx_temp=TemporalTkeField[8].getDataPtr();
    double *dtau_xydy_temp=TemporalTkeField[9].getDataPtr();
    double *dtau_yydy_temp=TemporalTkeField[10].getDataPtr();
    double *dtau_zydy_temp=TemporalTkeField[11].getDataPtr();
    double *dtau_xzdz_temp=TemporalTkeField[12].getDataPtr();
    double *dtau_yzdz_temp=TemporalTkeField[13].getDataPtr();
    double *dtau_zzdz_temp=TemporalTkeField[14].getDataPtr();
    
    //dtau_xxdx
    solveGrad(TemporalDerivedField[16], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_x(TemporalTkeField[6],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    //dtau_yxdx
    solveGrad(TemporalDerivedField[17], myOP, myBf[1], myMPI,myBC,1,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_x(TemporalTkeField[7],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    //dtau_zxdx
    solveGrad(TemporalDerivedField[18], myOP, myBf[2], myMPI,myBC,2,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_x(TemporalTkeField[8],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    //dtau_xydy
    solveGrad(TemporalDerivedField[17], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_y(TemporalTkeField[9],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    //dtau_yydy
    solveGrad(TemporalDerivedField[19], myOP, myBf[1], myMPI,myBC,1,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_y(TemporalTkeField[10],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    //dtau_zydy
    solveGrad(TemporalDerivedField[20], myOP, myBf[2], myMPI,myBC,2,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_y(TemporalTkeField[11],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    //dtau_xzdz
    solveGrad(TemporalDerivedField[18], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_z(TemporalTkeField[12],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    //dtau_yzdz
    solveGrad(TemporalDerivedField[20], myOP, myBf[1], myMPI,myBC,1,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_z(TemporalTkeField[13],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    //dtau_zzdz
    solveGrad(TemporalDerivedField[21], myOP, myBf[2], myMPI,myBC,2,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_z(TemporalTkeField[14],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    double *utau_xx_temp=TemporalTkeField[15].getDataPtr();
    double *vtau_yx_temp=TemporalTkeField[16].getDataPtr();
    double *wtau_zx_temp=TemporalTkeField[17].getDataPtr();
    double *uktau_kx_temp=TemporalTkeField[18].getDataPtr();
    double *duktau_kxdx_temp=AveragedTkeField[19].getDataPtr();
    
    double *utau_xy_temp=TemporalTkeField[20].getDataPtr();
    double *vtau_yy_temp=TemporalTkeField[21].getDataPtr();
    double *wtau_zy_temp=TemporalTkeField[22].getDataPtr();
    double *uktau_ky_temp=TemporalTkeField[23].getDataPtr();
    double *duktau_kydy_temp=AveragedTkeField[24].getDataPtr();
    
    double *utau_xz_temp=TemporalTkeField[25].getDataPtr();
    double *vtau_yz_temp=TemporalTkeField[26].getDataPtr();
    double *wtau_zz_temp=TemporalTkeField[27].getDataPtr();
    double *uktau_kz_temp=TemporalTkeField[28].getDataPtr();
    double *duktau_kzdz_temp=AveragedTkeField[29].getDataPtr();
    
    double *tau_xxdudx_temp=TemporalTkeField[30].getDataPtr();
    double *tau_yxdvdx_temp=TemporalTkeField[31].getDataPtr();
    double *tau_zxdwdx_temp=TemporalTkeField[32].getDataPtr();
    double *tau_xydudy_temp=TemporalTkeField[33].getDataPtr();
    double *tau_yydvdy_temp=TemporalTkeField[34].getDataPtr();
    double *tau_zydwdy_temp=TemporalTkeField[35].getDataPtr();
    double *tau_xzdudz_temp=TemporalTkeField[36].getDataPtr();
    double *tau_yzdvdz_temp=TemporalTkeField[37].getDataPtr();
    double *tau_zzdwdz_temp=TemporalTkeField[38].getDataPtr();
    
    double *rhouuu_temp=TemporalTkeField[39].getDataPtr();
    double *rhouuv_temp=TemporalTkeField[40].getDataPtr();
    double *rhouuw_temp=TemporalTkeField[41].getDataPtr();
    double *rhouvv_temp=TemporalTkeField[42].getDataPtr();
    double *rhouvw_temp=TemporalTkeField[43].getDataPtr();
    double *rhouww_temp=TemporalTkeField[44].getDataPtr();
    double *rhovvv_temp=TemporalTkeField[45].getDataPtr();
    double *rhovvw_temp=TemporalTkeField[46].getDataPtr();
    double *rhovww_temp=TemporalTkeField[47].getDataPtr();
    double *rhowww_temp=TemporalTkeField[48].getDataPtr();
    
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                
                up_temp[idx_xi]=u_temp[idx_xi]*p_temp[idx_xi];
                vp_temp[idx_xi]=v_temp[idx_xi]*p_temp[idx_xi];
                wp_temp[idx_xi]=w_temp[idx_xi]*p_temp[idx_xi];
                pdudx_temp[idx_xi]=p_temp[idx_xi]*ux_temp[idx_xi];
                pdvdy_temp[idx_xi]=p_temp[idx_xi]*vy_temp[idx_xi];
                pdwdz_temp[idx_xi]=p_temp[idx_xi]*wz_temp[idx_xi];
                
                utau_xx_temp[idx_xi]=u_temp[idx_xi]*tau_xx_temp[idx_xi];
                vtau_yx_temp[idx_xi]=v_temp[idx_xi]*tau_xy_temp[idx_xi];
                wtau_zx_temp[idx_xi]=w_temp[idx_xi]*tau_xz_temp[idx_xi];
                uktau_kx_temp[idx_xi]=utau_xx_temp[idx_xi]+vtau_yx_temp[idx_xi]+wtau_zx_temp[idx_xi];
                
                utau_xy_temp[idx_xi]=u_temp[idx_xi]*tau_xy_temp[idx_xi];
                vtau_yy_temp[idx_xi]=v_temp[idx_xi]*tau_yy_temp[idx_xi];
                wtau_zy_temp[idx_xi]=w_temp[idx_xi]*tau_yz_temp[idx_xi];
                uktau_ky_temp[idx_xi]=utau_xy_temp[idx_xi]+vtau_yy_temp[idx_xi]+wtau_zy_temp[idx_xi];
                
                utau_xz_temp[idx_xi]=u_temp[idx_xi]*tau_xz_temp[idx_xi];
                vtau_yz_temp[idx_xi]=v_temp[idx_xi]*tau_yz_temp[idx_xi];
                wtau_zz_temp[idx_xi]=w_temp[idx_xi]*tau_zz_temp[idx_xi];
                uktau_kz_temp[idx_xi]=utau_xz_temp[idx_xi]+vtau_yz_temp[idx_xi]+wtau_zz_temp[idx_xi];
                
                tau_xxdudx_temp[idx_xi]=tau_xx_temp[idx_xi]*ux_temp[idx_xi];
                tau_yxdvdx_temp[idx_xi]=tau_xy_temp[idx_xi]*vx_temp[idx_xi];
                tau_zxdwdx_temp[idx_xi]=tau_xz_temp[idx_xi]*wx_temp[idx_xi];
                tau_xydudy_temp[idx_xi]=tau_xy_temp[idx_xi]*ux_temp[idx_xi];
                tau_yydvdy_temp[idx_xi]=tau_yy_temp[idx_xi]*vx_temp[idx_xi];
                tau_zydwdy_temp[idx_xi]=tau_yz_temp[idx_xi]*wx_temp[idx_xi];
                tau_xzdudz_temp[idx_xi]=tau_xz_temp[idx_xi]*uz_temp[idx_xi];
                tau_yzdvdz_temp[idx_xi]=tau_yz_temp[idx_xi]*vz_temp[idx_xi];
                tau_zzdwdz_temp[idx_xi]=tau_zz_temp[idx_xi]*wz_temp[idx_xi];
                
                rhouuu_temp[idx_xi]=rho_temp[idx_xi]*u_temp[idx_xi]*u_temp[idx_xi]*u_temp[idx_xi];
                rhouuv_temp[idx_xi]=rho_temp[idx_xi]*u_temp[idx_xi]*u_temp[idx_xi]*v_temp[idx_xi];
                rhouuw_temp[idx_xi]=rho_temp[idx_xi]*u_temp[idx_xi]*w_temp[idx_xi];
                rhouvv_temp[idx_xi]=rho_temp[idx_xi]*u_temp[idx_xi]*v_temp[idx_xi]*v_temp[idx_xi];
                rhouvw_temp[idx_xi]=rho_temp[idx_xi]*u_temp[idx_xi]*v_temp[idx_xi]*w_temp[idx_xi];
                rhouww_temp[idx_xi]=rho_temp[idx_xi]*u_temp[idx_xi]*w_temp[idx_xi]*w_temp[idx_xi];
                rhovvv_temp[idx_xi]=rho_temp[idx_xi]*v_temp[idx_xi]*v_temp[idx_xi]*v_temp[idx_xi];
                rhovvw_temp[idx_xi]=rho_temp[idx_xi]*v_temp[idx_xi]*v_temp[idx_xi]*w_temp[idx_xi];
                rhovww_temp[idx_xi]=rho_temp[idx_xi]*v_temp[idx_xi]*w_temp[idx_xi]*w_temp[idx_xi];
                rhowww_temp[idx_xi]=rho_temp[idx_xi]*w_temp[idx_xi]*w_temp[idx_xi]*w_temp[idx_xi];
            }
        }
    }
    
    //duktau_kxdx
    solveGrad(TemporalTkeField[18], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_x(TemporalTkeField[19],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    //duktau_kydy
    solveGrad(TemporalTkeField[23], myOP, myBf[1], myMPI,myBC,1,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_y(TemporalTkeField[24],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    //duktau_kzdz
    solveGrad(TemporalTkeField[28], myOP, myBf[2], myMPI,myBC,2,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_z(TemporalTkeField[29],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    for(int n=0;n<VarNameList_TKE.size();n++)
    {
        double *f_aver=AveragedTkeField[n].getDataPtr();
        double *f_temp=TemporalTkeField[n].getDataPtr();
        
        for (int k=0; k<nz; k++)
        {
            for (int j=0; j<ny; j++)
            {
                for (int i=0; i<nx; i++)
                {
                    int idx_xi=nx*ny*k+nx*j+i;
                    
                    f_aver[idx_xi]=(f_aver[idx_xi]*averStep+f_temp[idx_xi])/(averStep+1.0);
                }
            }
        }
    }
    
    
    double *rho_aver=AveragedPrimaryField[0].getDataPtr();
    double *u_aver=AveragedPrimaryField[1].getDataPtr();
    double *v_aver=AveragedPrimaryField[2].getDataPtr();
    double *w_aver=AveragedPrimaryField[3].getDataPtr();
    double *p_aver=AveragedPrimaryField[4].getDataPtr();
    
    double* ux_aver=AveragedGradsField[0].getDataPtr();
    double* uy_aver=AveragedGradsField[1].getDataPtr();
    double* uz_aver=AveragedGradsField[2].getDataPtr();
    double* vx_aver=AveragedGradsField[3].getDataPtr();
    double* vy_aver=AveragedGradsField[4].getDataPtr();
    double* vz_aver=AveragedGradsField[5].getDataPtr();
    double* wx_aver=AveragedGradsField[6].getDataPtr();
    double* wy_aver=AveragedGradsField[7].getDataPtr();
    double* wz_aver=AveragedGradsField[8].getDataPtr();
    double *px_aver=AveragedGradsField[9].getDataPtr();
    double *py_aver=AveragedGradsField[10].getDataPtr();
    double *pz_aver=AveragedGradsField[11].getDataPtr();
    
    
    double *rhou_aver=AveragedDerivedField[0].getDataPtr();
    double *rhov_aver=AveragedDerivedField[1].getDataPtr();
    double *rhow_aver=AveragedDerivedField[2].getDataPtr();
    
    double *rhouu_aver=AveragedDerivedField[10].getDataPtr();
    double *rhouv_aver=AveragedDerivedField[11].getDataPtr();
    double *rhouw_aver=AveragedDerivedField[12].getDataPtr();
    double *rhovv_aver=AveragedDerivedField[13].getDataPtr();
    double *rhovw_aver=AveragedDerivedField[14].getDataPtr();
    double *rhoww_aver=AveragedDerivedField[15].getDataPtr();
    double *tau_xx_aver=AveragedDerivedField[16].getDataPtr();
    double *tau_xy_aver=AveragedDerivedField[17].getDataPtr();
    double *tau_xz_aver=AveragedDerivedField[18].getDataPtr();
    double *tau_yy_aver=AveragedDerivedField[19].getDataPtr();
    double *tau_yz_aver=AveragedDerivedField[20].getDataPtr();
    double *tau_zz_aver=AveragedDerivedField[21].getDataPtr();
    
    double *up_aver=AveragedTkeField[0].getDataPtr();
    double *vp_aver=AveragedTkeField[1].getDataPtr();
    double *wp_aver=AveragedTkeField[2].getDataPtr();
    double *pdudx_aver=AveragedTkeField[3].getDataPtr();
    double *pdvdy_aver=AveragedTkeField[4].getDataPtr();
    double *pdwdz_aver=AveragedTkeField[5].getDataPtr();
    
    double *dtau_xxdx_aver=AveragedTkeField[6].getDataPtr();
    double *dtau_yxdx_aver=AveragedTkeField[7].getDataPtr();
    double *dtau_zxdx_aver=AveragedTkeField[8].getDataPtr();
    double *dtau_xydy_aver=AveragedTkeField[9].getDataPtr();
    double *dtau_yydy_aver=AveragedTkeField[10].getDataPtr();
    double *dtau_zydy_aver=AveragedTkeField[11].getDataPtr();
    double *dtau_xzdz_aver=AveragedTkeField[12].getDataPtr();
    double *dtau_yzdz_aver=AveragedTkeField[13].getDataPtr();
    double *dtau_zzdz_aver=AveragedTkeField[14].getDataPtr();
    
    double *utau_xx_aver=AveragedTkeField[15].getDataPtr();
    double *vtau_yx_aver=AveragedTkeField[16].getDataPtr();
    double *wtau_zx_aver=AveragedTkeField[17].getDataPtr();
    double *uktau_kx_aver=AveragedTkeField[18].getDataPtr();
    double *duktau_kxdx_aver=AveragedTkeField[19].getDataPtr();
    
    double *utau_xy_aver=AveragedTkeField[20].getDataPtr();
    double *vtau_yy_aver=AveragedTkeField[21].getDataPtr();
    double *wtau_zy_aver=AveragedTkeField[22].getDataPtr();
    double *uktau_ky_aver=AveragedTkeField[23].getDataPtr();
    double *duktau_kydy_aver=AveragedTkeField[24].getDataPtr();
    
    double *utau_xz_aver=AveragedTkeField[25].getDataPtr();
    double *vtau_yz_aver=AveragedTkeField[26].getDataPtr();
    double *wtau_zz_aver=AveragedTkeField[27].getDataPtr();
    double *uktau_kz_aver=AveragedTkeField[28].getDataPtr();
    double *duktau_kzdz_aver=AveragedTkeField[29].getDataPtr();
    
    double *tau_xxdudx_aver=AveragedTkeField[30].getDataPtr();
    double *tau_yxdvdx_aver=AveragedTkeField[31].getDataPtr();
    double *tau_zxdwdx_aver=AveragedTkeField[32].getDataPtr();
    double *tau_xydudy_aver=AveragedTkeField[33].getDataPtr();
    double *tau_yydvdy_aver=AveragedTkeField[34].getDataPtr();
    double *tau_zydwdy_aver=AveragedTkeField[35].getDataPtr();
    double *tau_xzdudz_aver=AveragedTkeField[36].getDataPtr();
    double *tau_yzdvdz_aver=AveragedTkeField[37].getDataPtr();
    double *tau_zzdwdz_aver=AveragedTkeField[38].getDataPtr();
    
    double *tke=TKEbudgetField[0].getDataPtr();
    double *production=TKEbudgetField[1].getDataPtr();
    double *transportation=TKEbudgetField[2].getDataPtr();
    double *pressureDiffusion=TKEbudgetField[3].getDataPtr();
    double *p1=TKEbudgetField[4].getDataPtr();
    double *p2=TKEbudgetField[5].getDataPtr();
    double *p3=TKEbudgetField[6].getDataPtr();
    double *viscousDiffusion=TKEbudgetField[7].getDataPtr();
    double *v1=TKEbudgetField[8].getDataPtr();
    double *v2=TKEbudgetField[9].getDataPtr();
    double *v3=TKEbudgetField[10].getDataPtr();
    double *convection=TKEbudgetField[11].getDataPtr();
    
    double *drhoudx_aver=TKEbudgetField[12].getDataPtr();
    double *drhoudy_aver=TKEbudgetField[13].getDataPtr();
    double *drhoudz_aver=TKEbudgetField[14].getDataPtr();
    double *drhovdx_aver=TKEbudgetField[15].getDataPtr();
    double *drhovdy_aver=TKEbudgetField[16].getDataPtr();
    double *drhovdz_aver=TKEbudgetField[17].getDataPtr();
    double *drhowdx_aver=TKEbudgetField[18].getDataPtr();
    double *drhowdy_aver=TKEbudgetField[19].getDataPtr();
    double *drhowdz_aver=TKEbudgetField[20].getDataPtr();
    
    //drhou
    solveGrad(AveragedDerivedField[0], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(TKEbudgetField[12],TKEbudgetField[13],TKEbudgetField[14],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    //drhov
    solveGrad(AveragedDerivedField[1], myOP, myBf[1], myMPI,myBC,1,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(TKEbudgetField[15],TKEbudgetField[16],TKEbudgetField[17],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    //drhow
    solveGrad(AveragedDerivedField[2], myOP, myBf[2], myMPI,myBC,2,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_xyz(TKEbudgetField[18],TKEbudgetField[19],TKEbudgetField[20],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    double *u_p_aver=TKEbudgetField[21].getDataPtr();
    double *v_p_aver=TKEbudgetField[22].getDataPtr();
    double *w_p_aver=TKEbudgetField[23].getDataPtr();
    
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                u_p_aver[idx_xi]=(up_aver[idx_xi]-u_aver[idx_xi]*p_aver[idx_xi]);
                v_p_aver[idx_xi]=(vp_aver[idx_xi]-v_aver[idx_xi]*p_aver[idx_xi]);
                w_p_aver[idx_xi]=(wp_aver[idx_xi]-w_aver[idx_xi]*p_aver[idx_xi]);
            }
        }
    }
    //dupdx
    solveGrad(TKEbudgetField[21], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_x(TKEbudgetField[24],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    //dvpdy
    solveGrad(TKEbudgetField[22], myOP, myBf[1], myMPI,myBC,1,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_y(TKEbudgetField[25],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    //dwpdz
    solveGrad(TKEbudgetField[23], myOP, myBf[2], myMPI,myBC,2,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_z(TKEbudgetField[26],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    double *uk_taukx_aver=TKEbudgetField[27].getDataPtr();
    double *uk_tauky_aver=TKEbudgetField[28].getDataPtr();
    double *uk_taukz_aver=TKEbudgetField[29].getDataPtr();
    
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                
                uk_taukx_aver[idx_xi]=u_aver[idx_xi]*tau_xx_aver[idx_xi]
                +v_aver[idx_xi]*tau_xy_aver[idx_xi]
                +w_aver[idx_xi]*tau_xz_aver[idx_xi];
                
                uk_tauky_aver[idx_xi]=u_aver[idx_xi]*tau_xy_aver[idx_xi]
                +v_aver[idx_xi]*tau_yz_aver[idx_xi]
                +w_aver[idx_xi]*tau_yz_aver[idx_xi];
                
                uk_taukz_aver[idx_xi]=u_aver[idx_xi]*tau_xz_aver[idx_xi]
                +v_aver[idx_xi]*tau_yz_aver[idx_xi]
                +w_aver[idx_xi]*tau_zz_aver[idx_xi];
                
            }
        }
    }
    
    //dupdx
    solveGrad(TKEbudgetField[27], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_x(TKEbudgetField[30],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    //dvpdy
    solveGrad(TKEbudgetField[28], myOP, myBf[1], myMPI,myBC,1,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_y(TKEbudgetField[31],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    //dwpdz
    solveGrad(TKEbudgetField[29], myOP, myBf[2], myMPI,myBC,2,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_z(TKEbudgetField[32],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    
    double *dupdx=TKEbudgetField[24].getDataPtr();
    double *dvpdy=TKEbudgetField[25].getDataPtr();
    double *dwpdz=TKEbudgetField[26].getDataPtr();
    
    double *duktau_kxdx=TKEbudgetField[30].getDataPtr();
    double *duktau_kydy=TKEbudgetField[31].getDataPtr();
    double *duktau_kzdz=TKEbudgetField[32].getDataPtr();
    
    
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                
                tke[idx_xi]=0.5*(rhouu_aver[idx_xi]/rho_aver[idx_xi]-std::pow(rhou_aver[idx_xi]/rho_aver[idx_xi],2.0))
                +0.5*(rhovv_aver[idx_xi]/rho_aver[idx_xi]-std::pow(rhov_aver[idx_xi]/rho_aver[idx_xi],2.0))
                +0.5*(rhoww_aver[idx_xi]/rho_aver[idx_xi]-std::pow(rhow_aver[idx_xi]/rho_aver[idx_xi],2.0));
                
                production[idx_xi]=(rhouu_aver[idx_xi]-rhou_aver[idx_xi]*rhou_aver[idx_xi]/rho_aver[idx_xi])*drhoudx_aver[idx_xi]
                +(rhouv_aver[idx_xi]-rhou_aver[idx_xi]*rhov_aver[idx_xi]/rho_aver[idx_xi])*drhoudy_aver[idx_xi]
                +(rhouw_aver[idx_xi]-rhou_aver[idx_xi]*rhow_aver[idx_xi]/rho_aver[idx_xi])*drhoudz_aver[idx_xi]
                +(rhouv_aver[idx_xi]-rhou_aver[idx_xi]*rhov_aver[idx_xi]/rho_aver[idx_xi])*drhovdx_aver[idx_xi]
                +(rhovv_aver[idx_xi]-rhov_aver[idx_xi]*rhov_aver[idx_xi]/rho_aver[idx_xi])*drhovdy_aver[idx_xi]
                +(rhovw_aver[idx_xi]-rhov_aver[idx_xi]*rhow_aver[idx_xi]/rho_aver[idx_xi])*drhovdz_aver[idx_xi]
                +(rhouw_aver[idx_xi]-rhou_aver[idx_xi]*rhow_aver[idx_xi]/rho_aver[idx_xi])*drhowdx_aver[idx_xi]
                +(rhovw_aver[idx_xi]-rhov_aver[idx_xi]*rhow_aver[idx_xi]/rho_aver[idx_xi])*drhowdy_aver[idx_xi]
                +(rhoww_aver[idx_xi]-rhow_aver[idx_xi]*rhow_aver[idx_xi]/rho_aver[idx_xi])*drhowdz_aver[idx_xi];
                
                p1[idx_xi]=-(u_aver[idx_xi]-rhou_aver[idx_xi]/rho_aver[idx_xi])*px_aver[idx_xi]
                -(v_aver[idx_xi]-rhov_aver[idx_xi]/rho_aver[idx_xi])*py_aver[idx_xi]
                -(w_aver[idx_xi]-rhow_aver[idx_xi]/rho_aver[idx_xi])*pz_aver[idx_xi];
                
                p2[idx_xi]=-dupdx[idx_xi]-dvpdy[idx_xi]-dwpdz[idx_xi];
                
                p3[idx_xi]=(pdudx_aver[idx_xi]-ux_aver[idx_xi]*p_aver[idx_xi])
                +(pdvdy_aver[idx_xi]-vy_aver[idx_xi]*p_aver[idx_xi])
                +(pdwdz_aver[idx_xi]-wz_aver[idx_xi]*p_aver[idx_xi]);
                
                pressureDiffusion[idx_xi]=p1[idx_xi]+p2[idx_xi]+p3[idx_xi];
                
                v1[idx_xi]=(u_aver[idx_xi]-rhou_aver[idx_xi]/rho_aver[idx_xi])*dtau_xxdx_aver[idx_xi]
                +(v_aver[idx_xi]-rhov_aver[idx_xi]/rho_aver[idx_xi])*dtau_yxdx_aver[idx_xi]
                +(w_aver[idx_xi]-rhow_aver[idx_xi]/rho_aver[idx_xi])*dtau_zxdx_aver[idx_xi]
                +(u_aver[idx_xi]-rhou_aver[idx_xi]/rho_aver[idx_xi])*dtau_xydy_aver[idx_xi]
                +(v_aver[idx_xi]-rhov_aver[idx_xi]/rho_aver[idx_xi])*dtau_yydy_aver[idx_xi]
                +(w_aver[idx_xi]-rhow_aver[idx_xi]/rho_aver[idx_xi])*dtau_zydy_aver[idx_xi]
                +(u_aver[idx_xi]-rhou_aver[idx_xi]/rho_aver[idx_xi])*dtau_xzdz_aver[idx_xi]
                +(v_aver[idx_xi]-rhov_aver[idx_xi]/rho_aver[idx_xi])*dtau_yzdz_aver[idx_xi]
                +(w_aver[idx_xi]-rhow_aver[idx_xi]/rho_aver[idx_xi])*dtau_zzdz_aver[idx_xi];
                
                v2[idx_xi]=(duktau_kxdx_aver-duktau_kxdx)
                +(duktau_kydy_aver-duktau_kydy)
                +(duktau_kxdx_aver-duktau_kzdz);
                
                v3[idx_xi]=-(tau_xxdudx_aver[idx_xi]-ux_aver[idx_xi]*tau_xx_aver[idx_xi])
                -(tau_yxdvdx_aver[idx_xi]-vx_aver[idx_xi]*tau_xy_aver[idx_xi])
                -(tau_zxdwdx_aver[idx_xi]-wx_aver[idx_xi]*tau_xz_aver[idx_xi])
                -(tau_xydudy_aver[idx_xi]-uy_aver[idx_xi]*tau_xy_aver[idx_xi])
                -(tau_yydvdy_aver[idx_xi]-vy_aver[idx_xi]*tau_yy_aver[idx_xi])
                -(tau_zydwdy_aver[idx_xi]-wy_aver[idx_xi]*tau_yz_aver[idx_xi])
                -(tau_xzdudz_aver[idx_xi]-uz_aver[idx_xi]*tau_xz_aver[idx_xi])
                -(tau_yzdvdz_aver[idx_xi]-vz_aver[idx_xi]*tau_yz_aver[idx_xi])
                -(tau_zzdwdz_aver[idx_xi]-wz_aver[idx_xi]*tau_zz_aver[idx_xi]);
                
                viscousDiffusion[idx_xi]=v1[idx_xi]+v2[idx_xi]+v3[idx_xi];
                
            }
        }
    }
    
    double *rhouk=TKEbudgetField[33].getDataPtr();
    double *rhovk=TKEbudgetField[34].getDataPtr();
    double *rhowk=TKEbudgetField[35].getDataPtr();
    
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                
                rhouk[idx_xi]=rhou_aver[idx_xi]*tke[idx_xi];
                rhovk[idx_xi]=rhov_aver[idx_xi]*tke[idx_xi];
                rhowk[idx_xi]=rhow_aver[idx_xi]*tke[idx_xi];
                
            }
        }
    }
    
    //drhoukdx
    solveGrad(TKEbudgetField[33], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_x(TKEbudgetField[36],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    //drhovkdy
    solveGrad(TKEbudgetField[34], myOP, myBf[1], myMPI,myBC,1,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_y(TKEbudgetField[37],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    //drhowkdz
    solveGrad(TKEbudgetField[35], myOP, myBf[2], myMPI,myBC,2,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_z(TKEbudgetField[38],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    double *drhoukdx=TKEbudgetField[36].getDataPtr();
    double *drhovkdy=TKEbudgetField[37].getDataPtr();
    double *drhowkdz=TKEbudgetField[38].getDataPtr();
    
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                convection[idx_xi]=drhoukdx[idx_xi]+drhovkdy[idx_xi]+drhowkdz[idx_xi];
            }
        }
    }
    
    double *u_tke=TKEbudgetField[39].getDataPtr();
    double *v_tke=TKEbudgetField[40].getDataPtr();
    double *w_tke=TKEbudgetField[41].getDataPtr();
    
    double *rhouuu_aver=AveragedTkeField[39].getDataPtr();
    double *rhouuv_aver=AveragedTkeField[40].getDataPtr();
    double *rhouuw_aver=AveragedTkeField[41].getDataPtr();
    double *rhouvv_aver=AveragedTkeField[42].getDataPtr();
    double *rhouvw_aver=AveragedTkeField[43].getDataPtr();
    double *rhouww_aver=AveragedTkeField[44].getDataPtr();
    double *rhovvv_aver=AveragedTkeField[45].getDataPtr();
    double *rhovvw_aver=AveragedTkeField[46].getDataPtr();
    double *rhovww_aver=AveragedTkeField[47].getDataPtr();
    double *rhowww_aver=AveragedTkeField[48].getDataPtr();
    
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                
                u_tke[idx_xi]=0.5*rhouuu_aver[idx_xi]-rhouu_aver[idx_xi]*rhou_aver[idx_xi]/rho_aver[idx_xi]+0.5*rhou_aver[idx_xi]*rhou_aver[idx_xi]*rhou_aver[idx_xi]/rho_aver[idx_xi]/rho_aver[idx_xi]
                +0.5*rhouvv_aver[idx_xi]-rhouv_aver[idx_xi]*rhov_aver[idx_xi]/rho_aver[idx_xi]+0.5*rhou_aver[idx_xi]*rhov_aver[idx_xi]*rhov_aver[idx_xi]/rho_aver[idx_xi]/rho_aver[idx_xi]
                +0.5*rhouww_aver[idx_xi]-rhouw_aver[idx_xi]*rhow_aver[idx_xi]/rho_aver[idx_xi]+0.5*rhou_aver[idx_xi]*rhow_aver[idx_xi]*rhow_aver[idx_xi]/rho_aver[idx_xi]/rho_aver[idx_xi]
                -rhou_aver[idx_xi]/rho_aver[idx_xi]*0.5*(rhouu_aver[idx_xi]-rhou_aver[idx_xi]*rhou_aver[idx_xi]/rho_aver[idx_xi]
                                                         +rhovv_aver[idx_xi]-rhov_aver[idx_xi]*rhov_aver[idx_xi]/rho_aver[idx_xi]
                                                         +rhoww_aver[idx_xi]-rhow_aver[idx_xi]*rhow_aver[idx_xi]/rho_aver[idx_xi]);
                
                v_tke[idx_xi]=0.5*rhouuv_aver[idx_xi]-rhouv_aver[idx_xi]*rhou_aver[idx_xi]/rho_aver[idx_xi]+0.5*rhov_aver[idx_xi]*rhou_aver[idx_xi]*rhou_aver[idx_xi]/rho_aver[idx_xi]/rho_aver[idx_xi]
                +0.5*rhovvv_aver[idx_xi]-rhovv_aver[idx_xi]*rhov_aver[idx_xi]/rho_aver[idx_xi]+0.5*rhov_aver[idx_xi]*rhov_aver[idx_xi]*rhov_aver[idx_xi]/rho_aver[idx_xi]/rho_aver[idx_xi]
                +0.5*rhovww_aver[idx_xi]-rhovw_aver[idx_xi]*rhow_aver[idx_xi]/rho_aver[idx_xi]+0.5*rhov_aver[idx_xi]*rhow_aver[idx_xi]*rhow_aver[idx_xi]/rho_aver[idx_xi]/rho_aver[idx_xi]
                -rhov_aver[idx_xi]/rho_aver[idx_xi]*0.5*(rhouu_aver[idx_xi]-rhou_aver[idx_xi]*rhou_aver[idx_xi]/rho_aver[idx_xi]
                                                         +rhovv_aver[idx_xi]-rhov_aver[idx_xi]*rhov_aver[idx_xi]/rho_aver[idx_xi]
                                                         +rhoww_aver[idx_xi]-rhow_aver[idx_xi]*rhow_aver[idx_xi]/rho_aver[idx_xi]);
                
                w_tke[idx_xi]=0.5*rhouuw_aver[idx_xi]-rhouw_aver[idx_xi]*rhou_aver[idx_xi]/rho_aver[idx_xi]+0.5*rhow_aver[idx_xi]*rhou_aver[idx_xi]*rhou_aver[idx_xi]/rho_aver[idx_xi]/rho_aver[idx_xi]
                +0.5*rhovvw_aver[idx_xi]-rhovw_aver[idx_xi]*rhov_aver[idx_xi]/rho_aver[idx_xi]+0.5*rhow_aver[idx_xi]*rhov_aver[idx_xi]*rhov_aver[idx_xi]/rho_aver[idx_xi]/rho_aver[idx_xi]
                +0.5*rhowww_aver[idx_xi]-rhoww_aver[idx_xi]*rhow_aver[idx_xi]/rho_aver[idx_xi]+0.5*rhow_aver[idx_xi]*rhow_aver[idx_xi]*rhow_aver[idx_xi]/rho_aver[idx_xi]/rho_aver[idx_xi]
                -rhow_aver[idx_xi]/rho_aver[idx_xi]*0.5*(rhouu_aver[idx_xi]-rhou_aver[idx_xi]*rhou_aver[idx_xi]/rho_aver[idx_xi]
                                                         +rhovv_aver[idx_xi]-rhov_aver[idx_xi]*rhov_aver[idx_xi]/rho_aver[idx_xi]
                                                         +rhoww_aver[idx_xi]-rhow_aver[idx_xi]*rhow_aver[idx_xi]/rho_aver[idx_xi]);
                
                
            }
        }
    }
    
    
    //drhoukdx
    solveGrad(TKEbudgetField[39], myOP, myBf[0], myMPI,myBC,0,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_x(TKEbudgetField[42],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    //drhovkdy
    solveGrad(TKEbudgetField[40], myOP, myBf[1], myMPI,myBC,1,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_y(TKEbudgetField[43],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    //drhowkdz
    solveGrad(TKEbudgetField[41], myOP, myBf[2], myMPI,myBC,2,dtempdxi,dtempdeta,dtempdzeta);
    solveGrad_z(TKEbudgetField[44],dtempdxi,dtempdeta,dtempdzeta,xi_xyz, eta_xyz, zeta_xyz);
    
    double *du_tkedx=TKEbudgetField[42].getDataPtr();
    double *dv_tkedy=TKEbudgetField[43].getDataPtr();
    double *dw_tkedz=TKEbudgetField[44].getDataPtr();
    
    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                transportation[idx_xi]=du_tkedx[idx_xi]+dv_tkedy[idx_xi]+dw_tkedz[idx_xi];
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
    solveGrad_df(temp_xi, myOP, myBf, myMPI, dfdxi, 0, fdID, typeL, typeR);
    
    typeL=myBC.getBCtype(2);
    typeR=myBC.getBCtype(3);
    solveGrad_df(temp_eta, myOP, myBf, myMPI, dfdeta, 1, fdID, typeL, typeR);
    
    typeL=myBC.getBCtype(4);
    typeR=myBC.getBCtype(5);
    solveGrad_df(temp_zeta, myOP, myBf, myMPI, dfdzeta, 2, fdID, typeL, typeR);
    
    
}

void nuc3d::postproc::solveGrad_x(Field Q_x,
                                  Field Q_xi,
                                  Field Q_eta,
                                  Field Q_zeta,
                                  VectorField &xi_xyz,
                                  VectorField &eta_xyz,
                                  VectorField &zeta_xyz)
{
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
    double *pf_xi=temp_xi.getDataPtr();
    double *pf_eta=temp_eta.getDataPtr();
    double *pf_zeta=temp_zeta.getDataPtr();
    
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

void nuc3d::postproc::writeField_binary(std::ofstream &myFile, nuc3d::Field &myField)
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


void nuc3d::postproc::readField_binary(std::ifstream &myFile, nuc3d::Field &myField)
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


