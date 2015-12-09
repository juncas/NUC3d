//
//  NaiverStokesData3d.hpp
//  NUC3d
//
//  Created by Jun Peng on 15/11/3.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#ifndef NaiverStokesData3d_hpp
#define NaiverStokesData3d_hpp
#include "euler3d.h"


namespace nuc3d
{
    class gradvector
    {
        friend class NaiverStokesData3d;
        Field dxi;
        Field deta;
        Field dzeta;
    public:
        gradvector(int,int,int);
        ~gradvector();
    private:
        Field& getdxi(){return dxi;};
        Field& getdeta(){return deta;};
        Field& getdzeta(){return dzeta;};
    };
    
    class NaiverStokesData3d : virtual public EulerData3D
    {
        friend class physicsModel;
    public:
        gradvector du;
        gradvector dv;
        gradvector dw;
        gradvector dT;
        
        Field miu;
        Field coeff;
        
        //viscous stresses
        VectorField tau;
        
        //viscous fluxes
        VectorField Flux_xi_vis;
        VectorField Flux_eta_vis;
        VectorField Flux_zeta_vis;
        
        //viscous derivatives
        VectorField dfvdxi;
        VectorField dgvdeta;
        VectorField dhvdzeta;
    public:
        NaiverStokesData3d(int,int,int,int);
        ~NaiverStokesData3d();
        
    public:
        
        virtual void solve(PDEData3d &,
                           fieldOperator3d &,
                           VectorBuffer &,
                           physicsModel &,
                           MPIComunicator3d_nonblocking &,
                           boundaryCondition &);
        
        
    protected:
        
        void solveVis(PDEData3d &,
                      fieldOperator3d &,
                      physicsModel &,
                      std::vector<bufferData> &,
                      MPIComunicator3d_nonblocking &,
                      boundaryCondition &);
        
    private:
        void setBoundaryGrad(PDEData3d &myPDE,
                             fieldOperator3d &myOP,
                             physicsModel &myModel,
                             std::vector<bufferData> &myBf,
                             MPIComunicator3d_nonblocking &myMPI,
                             boundaryCondition &myBC);
        
        void solveGrad();
        void solveGradXi();
        void solveGradEta();
        void solveGradZeta();
        void solveStress();
        
        void setBoundaryViscousFlux();
        void solveViscousFlux();
        
        void setDerivativesVis();
        void solveRHS(PDEData3d &);
        
    };
}
#endif /* NaiverStokesData3d_hpp */
