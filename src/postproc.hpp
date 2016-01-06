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
    class IOController;

    
    class postproc
    {
        int nx;
        int ny;
        int nz;
        int postSize;
        Field Q;
        Field omega0;
        Field omega1;
        Field omega2;
        /*
         Other statistics can be added
         */
        
        Field q_xi,q_eta,q_zeta;
        
        Field u_xi,u_eta,u_zeta;
        Field v_xi,v_eta,v_zeta;
        Field w_xi,w_eta,w_zeta;
        
        Field u_x,u_y,u_z;
        Field v_x,v_y,v_z;
        Field w_x,w_y,w_z;
        
        double enstrophy_glb,enstrophy;
        
    public:
        postproc(int,int,int);
        ~postproc();
        
        void solvePost(VectorField &prims,
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
        
        void solveGrad_xyz(VectorField &xi_xyz,
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
        void writeField(std::ofstream &myFile, Field &myField);

    };
}

#endif /* postproc_hpp */
