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
    class NaiverStokesData3d : virtual public EulerData3D
    {
        friend class physicsModel;
    public:
        VectorField du;
        VectorField dv;
        VectorField dw;
        VectorField dT;
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
        virtual VectorField& getDrivativeXi();
        virtual VectorField& getDrivativeEta();
        virtual VectorField& getDrivativeZeta();
        virtual void solveLocal();
        
        virtual void solve(PDEData3d &,
                           fieldOperator3d &,
                           bufferData &,
                           physicsModel &,
                           MPIComunicator3d_nonblocking &);
        
    protected:
        
        void solveVis(fieldOperator3d &,
                      physicsModel &,
                      bufferData &,
                      MPIComunicator3d_nonblocking &);
        void setViscousFluxes();
    private:
        virtual void solveRHS(PDEData3d &);
        
    };
}
#endif /* NaiverStokesData3d_hpp */
