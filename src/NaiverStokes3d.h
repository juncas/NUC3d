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
#include "gradvector.hpp"


namespace nuc3d
{    
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
        
        VectorField &getVisFlux_xi(){return Flux_xi_vis;};
        VectorField &getVisFlux_eta(){return Flux_eta_vis;};
        VectorField &getVisFlux_zeta(){return Flux_zeta_vis;};
        
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
       
        void solveGrads(PDEData3d &myPDE,
                        fieldOperator3d &myOP,
                        std::vector<bufferData> &myBf,
                        MPIComunicator3d_nonblocking &myMPI,
                        boundaryCondition &myBC);
        
        void solve_grad(Field &,
                        fieldOperator3d &,
                        bufferData &,
                        MPIComunicator3d_nonblocking &,
                        gradvector &,
                        int fdID,
                        boundaryCondition &myBC);
        
        void solveGradXi(Field &,
                         fieldOperator3d &,
                         bufferData &myBf,
                         MPIComunicator3d_nonblocking &,
                         Field &,
                         int fdID,int typeL,int typeR);
        
        void solveGradEta(Field &,
                          fieldOperator3d &,
                          bufferData &myBf,
                          MPIComunicator3d_nonblocking &,
                          Field &,
                          int fdID,int typeL,int typeR);
        
        void solveGradZeta(Field &,
                           fieldOperator3d &,
                           bufferData &myBf,
                           MPIComunicator3d_nonblocking &,
                           Field &,
                           int fdID,int typeL,int typeR);
        
        void solveViscousFlux(physicsModel &myPhyMod);
                
        void setDerivativesVis(fieldOperator3d &myOP,
                               std::vector<bufferData> &myBf,
                               MPIComunicator3d_nonblocking &myMPI,
                               boundaryCondition &myBC);
        
        void setDerivativeXi(fieldOperator3d &myOP,
                             std::vector<bufferData> &myBf,
                             MPIComunicator3d_nonblocking &myMPI,
                             boundaryCondition &myBC);
        
        void setDerivativeEta(fieldOperator3d &myOP,
                              std::vector<bufferData> &myBf,
                              MPIComunicator3d_nonblocking &myMPI,
                              boundaryCondition &myBC);
        
        void setDerivativeZeta(fieldOperator3d &myOP,
                               std::vector<bufferData> &myBf,
                               MPIComunicator3d_nonblocking &myMPI,
                               boundaryCondition &myBC);
        
        void solveRHS(PDEData3d &);
        
    };
}
#endif /* NaiverStokesData3d_hpp */
