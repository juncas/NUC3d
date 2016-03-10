//
//  postproc.hpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/4.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//

#ifndef postproc_hpp
#define postproc_hpp

#include <stdio.h>
#include "field.h"
#include "mpi.h"
#include "bufferData.hpp"


namespace nuc3d
{
    class fieldOperator3d;
    class bufferData;
    class MPIComunicator3d_nonblocking;
    class boundaryCondition;
    class physicsModel;
    class IOController;
    
    
    
    class postproc
    {
        int nx;
        int ny;
        int nz;
        int postSize;
        
        typedef void (postproc::*pSolveAveraged)(VectorField &prims,
                                                 VectorField &acous,
                                                 VectorField &xi_xyz,
                                                 VectorField &eta_xyz,
                                                 VectorField &zeta_xyz,
                                                 physicsModel &myPhys,
                                                 fieldOperator3d &myOP,
                                                 VectorBuffer &myBf,
                                                 MPIComunicator3d_nonblocking &myMPI,
                                                 boundaryCondition &myBC);
        
        const std::vector<std::string> VarNameList_primary={
            "rho","u","v","w","p","T","e","alpha","miu","coeff"       };
        
        const std::vector<std::string> VarNameList_grads={
            "dudx","dudy","dudz",
            "dvdx","dvdy","dvdz",
            "dwdx","dwdy","dwdz",
            "dpdx","dpdy","dpdz",
            "dTdx","dTdy","dTdz",
            "omega0","omega1","omega2",
            "omega0^2","omega1^2","omega2^2"
        };
        
        
        const std::vector<std::string> VarNameList_derived={
            "rhou","rhov","rhow","rhoe",
            "uu","uv","uw","vv","vw","ww",
            "rhouu","rhouv","rhouw","rhovv","rhovw","rhoww",
            "tau_xx","tau_xy","tau_xz","tau_yy","tau_yz","tau_zz",
            "rhotau_xx","rhotau_xy","rhotau_xz","rhotau_yy","rhotau_yz","rhotau_zz"
        };
        
        const std::vector<std::string> VarNameList_TKE={
            "up","vp","wp",
            "pdudx","pdvdy","pdwdz",
            "dtau_xxdx","dtau_yxdx","dtau_zxdx",
            "dtau_xydy","dtau_yydy","dtau_zydy",
            "dtau_xzdz","dtau_yzdz","dtau_zzdz",
            "utau_xx","vtau_yx","wtau_zx","uktau_kx","duktau_kxdx",
            "utau_xy","vtau_yy","wtau_zy","uktau_ky","duktau_kydy",
            "utau_xz","vtau_yz","wtau_zz","uktau_kz","duktau_kzdz",
            "tau_xxdudx","tau_yxdvdx","tau_zxdwdx",
            "tau_xydudy","tau_yydvdy","tau_zydwdy",
            "tau_xzdudz","tau_yzdvdz","tau_zzdwdz",
            "rhouuu","rhouuv","rhouuw",
            "rhovuv","rhouvw",
            "rhouww","rhovvv","rhovvw",
            "rhovww","rhowww",
        };
        
        const std::vector<std::string> VarNameList_TKEbudget={
            "tke",
            "production",
            "turbulent transportation",
            "pressure diffusion",
            "p1",
            "p2",
            "p3",
            "viscous diffusion",
            "v1",
            "v2",
            "v3",
            "convection",
            "drhoudx","drhoudy","drhoudz",
            "drhovdx","drhovdy","drhovdz",
            "drhowdx","drhowdy","drhowdz",
            "u'p'","v'p'","w'p'",
            "dupdx","dvpdy","dwpdz",
            "uk'tau_kx'","uk'tau_kx'","uk'tau_kx'",
            "duk'tau_kx'dx","duk'tau_kx'dx","duk'tau_kx'dx",
            "rhouk","rhovk","rhowk",
            "drhoukdx","drhovkdy","drhowkdz",
            "u_tke","v_tke","w_tke",
            "du_tkedx","dv_tkedy","dw_tkedz"
        };
        
        const std::vector<std::string> VarNameList_q={
            "dudx","dudy","dudz",
            "dvdx","dvdy","dvdz",
            "dwdx","dwdy","dwdz",
            "omega0","omega1","omega2",
            "enstrophy",
            "Q"
        };// other kinds of simutaneous variables can be added
        
        
        
        Field temp_xi,temp_eta,temp_zeta;
        Field dtempdxi,dtempdeta,dtempdzeta;
        
        VectorField TemporalPrimaryField;
        VectorField TemporalGradsField;
        VectorField TemporalDerivedField;
        VectorField TemporalTkeField;
        VectorField TemporalQField;
        
        VectorField AveragedPrimaryField;
        VectorField AveragedGradsField;
        VectorField AveragedDerivedField;
        VectorField AveragedTkeField;
        
        VectorField TKEbudgetField;
        
        double enstrophy_glb,enstrophy;
        double kinetic_glb,kinetic;
        
        std::vector<int> variableScalarInt;
        std::vector<double> variableScalarDouble;
        
        int averStep;
        
        pSolveAveraged myAverage[3]={
            &postproc::solveAveraged0,
            &postproc::solveAveraged1,
            &postproc::solveAveraged2
        };
        
    public:
        postproc(int,int,int);
        ~postproc();
        
        void initPost();
        
        void solvePost(VectorField &prims,
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
                       double time);
        
        void OutputPost(VectorField &prims,
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
                        double time);
    private:
        void solveGrad(Field &myField,
                       fieldOperator3d &myOP,
                       bufferData &myBf,
                       MPIComunicator3d_nonblocking &myMPI,
                       boundaryCondition &myBC,
                       int fdID,
                       Field &dfdxi,
                       Field &dfdeta,
                       Field &dfdzeta);
        
        void solveGrad_xyz(Field Q_x,
                           Field Q_y,
                           Field Q_z,
                           Field Q_xi,
                           Field Q_eta,
                           Field Q_zeta,
                           VectorField &xi_xyz,
                           VectorField &eta_xyz,
                           VectorField &zeta_xyz);
        
        void solveGrad_x(Field Q_x,
                         Field Q_xi,
                         Field Q_eta,
                         Field Q_zeta,
                         VectorField &xi_xyz,
                         VectorField &eta_xyz,
                         VectorField &zeta_xyz);
        
        void solveGrad_y(Field Q_y,
                         Field Q_xi,
                         Field Q_eta,
                         Field Q_zeta,
                         VectorField &xi_xyz,
                         VectorField &eta_xyz,
                         VectorField &zeta_xyz);
        
        void solveGrad_z(Field Q_z,
                         Field Q_xi,
                         Field Q_eta,
                         Field Q_zeta,
                         VectorField &xi_xyz,
                         VectorField &eta_xyz,
                         VectorField &zeta_xyz);
        
        void setField(const Field &f);
        
        void solveGrad_df(Field &myField,
                          fieldOperator3d &myOP,
                          bufferData &myBf,
                          MPIComunicator3d_nonblocking &myMPI,
                          Field &df,
                          int dir,
                          int fdID,
                          int typeL,
                          int typeR);
        
        void solveTemporal(VectorField &prims,
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
                           double time);
        
        void solveQ(VectorField &prims,
                    VectorField &xi_xyz,
                    VectorField &eta_xyz,
                    VectorField &zeta_xyz,
                    fieldOperator3d &myOP,
                    VectorBuffer &myBf,
                    MPIComunicator3d_nonblocking &myMPI,
                    boundaryCondition &myBC,
                    int istep,
                    double time);
        
        void solveAveraged(VectorField &prims,
                           VectorField &acous,
                           VectorField &xi_xyz,
                           VectorField &eta_xyz,
                           VectorField &zeta_xyz,
                           physicsModel &myPhys,
                           fieldOperator3d &myOP,
                           VectorBuffer &myBf,
                           MPIComunicator3d_nonblocking &myMPI,
                           boundaryCondition &myBC);
        
        void solveAveraged0(VectorField &prims,
                            VectorField &acous,
                            VectorField &xi_xyz,
                            VectorField &eta_xyz,
                            VectorField &zeta_xyz,
                            physicsModel &myPhys,
                            fieldOperator3d &myOP,
                            VectorBuffer &myBf,
                            MPIComunicator3d_nonblocking &myMPI,
                            boundaryCondition &myBC);
        
        void solveAveraged1(VectorField &prims,
                            VectorField &acous,
                            VectorField &xi_xyz,
                            VectorField &eta_xyz,
                            VectorField &zeta_xyz,
                            physicsModel &myPhys,
                            fieldOperator3d &myOP,
                            VectorBuffer &myBf,
                            MPIComunicator3d_nonblocking &myMPI,
                            boundaryCondition &myBC);
        
        void solveAveraged2(VectorField &prims,
                            VectorField &acous,
                            VectorField &xi_xyz,
                            VectorField &eta_xyz,
                            VectorField &zeta_xyz,
                            physicsModel &myPhys,
                            fieldOperator3d &myOP,
                            VectorBuffer &myBf,
                            MPIComunicator3d_nonblocking &myMPI,
                            boundaryCondition &myBC);
        
        void writeField(std::ofstream &myFile, Field &myField);
        void writeField_binary(std::ofstream &myFile, nuc3d::Field &myField);
        void readField_binary(std::ifstream &myFile, nuc3d::Field &myField);
        
        void OutputTemporal(VectorField &prims,
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
                            double time);
        
        void OutputAveraged_bin(VectorField &prims,
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
                            double time);
        
        void OutputAveraged_tec(VectorField &prims,
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
                            double time);

        
    };
}

#endif /* postproc_hpp */
