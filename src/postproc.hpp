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
        
        Field temp_xi,temp_eta,temp_zeta;
        Field dtempdxi,dtempdeta,dtempdzeta;
        
        VectorField variableField;
        std::vector<int> variableScalarInt;
        std::vector<double> variableScalarDouble;
        //
        //        Field u_x,u_y,u_z;
        //        Field v_x,v_y,v_z;
        //        Field w_x,w_y,w_z;
        //        Field p_x,p_y,p_z;
        //        Field miu,coeff;
        
        //
        //        Field Q;
        //
        //        Field omega0;
        //        Field omega1;
        //        Field omega2;
        //
        //        Field tau_xx,tau_xy,tau_xz;
        //        Field tau_yy,tau_yz;
        //        Field tau_zz;
        //
        //
        //
        //        Field rhou,rhov,rhow;
        //
        //        Field tau_xx_x,tau_xx_y,tau_xx_z;
        //        Field tau_xy_x,tau_xy_y,tau_xy_z;
        //        Field tau_xz_x,tau_xz_y,tau_xz_z;
        //        Field tau_yy_x,tau_yy_y,tau_yy_z;
        //        Field tau_yz_x,tau_yz_y,tau_yz_z;
        //        Field tau_zz_x,tau_zz_y,tau_zz_z;
        
        //        Field drhoudx,drhoudy,drhoudz;
        //        Field drhovdx,drhovdy,drhovdz;
        //        Field drhowdx,drhowdy,drhowdz;
        //
        //        /*turbulent statistic variables*/
        //        Field aver_rho,aver_u,aver_v,aver_w,aver_p,aver_T;
        //
        //
        //        //TKE budaget terms
        //        Field TKE;
        //
        //
        //        //Production
        //        Field P_tke;
        //        Field aver_rhou,aver_rhov,aver_rhow;
        //        Field aver_rhouu,aver_rhouv,aver_rhouw;
        //        Field aver_rhovv,aver_rhovw;
        //        Field aver_rhoww;
        //
        //        Field daver_rhoudx,daver_rhoudy,daver_rhoudz;
        //        Field daver_rhovdx,daver_rhovdy,daver_rhovdz;
        //        Field daver_rhowdx,daver_rhowdy,daver_rhowdz;
        //
        //        //Turbulent transportation
        //        Field T_tke;
        //        Field dTdx,dTdy,dTdz;
        //        Field aver_rhouuu,aver_rhouuv,aver_rhouuw;
        //        Field aver_rhovvu,aver_rhovvv,aver_rhovvw;
        //        Field aver_rhowwu,aver_rhowwv,aver_rhowww;
        //
        //        //Pressure dilatation
        //        Field PIE_tke;
        //        Field PIE_tke_0,PIE_tke_1,PIE_tke_2;
        //        Field aver_dpdx,aver_dpdy,aver_dpdz;
        //        Field aver_dpudx,aver_dpvdy,aver_dpwdz;
        //        Field aver_pu,aver_pv,aver_pw;
        //        Field aver_pdudx,aver_pdvdy,aver_pdwdz;
        //
        //        //mean stress work by fluctuation
        //        Field M_tke;
        //
        //        Field aver_tau_xx,aver_tau_xy,aver_tau_xz;
        //        Field aver_tau_yy,aver_tau_yz;
        //        Field aver_tau_zz;
        //
        //        Field aver_dtau_xxdx,aver_dtau_xydy,aver_dtau_xzdz;
        //        Field aver_dtau_yxdx,aver_dtau_yydy,aver_dtau_yzdz;
        //        Field aver_dtau_zxdx,aver_dtau_zydy,aver_dtau_zzdz;
        //
        //        //viscous diffusion
        //        Field D_tke;
        //        Field aver_utau_xx,aver_vtau_yx,aver_wtau_zx;
        //        Field aver_utau_xy,aver_vtau_yy,aver_wtau_zy;
        //        Field aver_utau_xz,aver_vtau_yz,aver_wtau_zz;
        //
        //        //viscous dissipation
        //        Field EPSL_tke;
        //
        //        Field aver_tau_xxdudx,aver_tau_xydudy,aver_tau_xzdudz;
        //        Field aver_tau_yxdvdx,aver_tau_yydvdy,aver_tau_yzdvdz;
        //        Field aver_tau_zxdwdx,aver_tau_zydwdy,aver_tau_zzdwdz;
        //
        //        Field aver_dudx,aver_dudy,aver_dudz;
        //        Field aver_dvdx,aver_dvdy,aver_dvdz;
        //        Field aver_dwdx,aver_dwdy,aver_dwdz;
        //
        //        Field EPSL_tke_s,EPSL_tke_d,EPSL_tke_i;
        //        Field aver_miuomega02,aver_miuomega12,aver_miuomega22;
        //        Field aver_omega0,aver_omega1,aver_omega2;
        //        Field aver_rhoomega0,aver_rhoomega1,aver_rhoomega2;
        
        
        /*
         Other statistics can be added
         */
        
        
        
        
        double enstrophy_glb,enstrophy;
        double kinetic_glb,kinetic;
        
        
        int averStep;
        
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
                       physicsModel myPhys,
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
        
        void solveQ(VectorField &prims,
                    VectorField &xi_xyz,
                    VectorField &eta_xyz,
                    VectorField &zeta_xyz,
                    fieldOperator3d &myOP,
                    VectorBuffer &myBf,
                    MPIComunicator3d_nonblocking &myMPI,
                    boundaryCondition &myBC);
        
        void solveTemporal(VectorField &prims,
                           VectorField &acous,
                           VectorField &xi_xyz,
                           VectorField &eta_xyz,
                           VectorField &zeta_xyz,
                           physicsModel myPhys,
                           fieldOperator3d &myOP,
                           VectorBuffer &myBf,
                           MPIComunicator3d_nonblocking &myMPI,
                           boundaryCondition &myBC);
        
        void solveAveraged(VectorField &prims,
                           VectorField &acous,
                           VectorField &xi_xyz,
                           VectorField &eta_xyz,
                           VectorField &zeta_xyz,
                           physicsModel myPhys,
                           fieldOperator3d &myOP,
                           VectorBuffer &myBf,
                           MPIComunicator3d_nonblocking &myMPI,
                           boundaryCondition &myBC);
        
        void solveTKE(VectorField &prims,
                      VectorField &acous,
                      VectorField &xi_xyz,
                      VectorField &eta_xyz,
                      VectorField &zeta_xyz,
                      physicsModel myPhys,
                      fieldOperator3d &myOP,
                      VectorBuffer &myBf,
                      MPIComunicator3d_nonblocking &myMPI,
                      boundaryCondition &myBC);
        
        void writeField(std::ofstream &myFile, Field &myField);
        
    };
}

#endif /* postproc_hpp */
