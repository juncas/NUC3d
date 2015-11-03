#ifndef Euler3d_h
#define Euler3d_h

#include <cstdlib>
#include <iostream>
#include "physicsModel.h"
#include "fieldOperator.h"
#include "MPICommunicator.h"

namespace nuc3d
{
    class EulerData3D;
    class EulerFlux;
    class PDEData3d;
    
    class EulerFlux
    {
        friend class physicsModel;
        friend class singleBlock;
        friend class EulerData3D;
        
        VectorField FluxL;
        VectorField FluxR;
        VectorField reconstFluxL;
        VectorField reconstFluxR;
        VectorField reconstFlux;
        double maxEigen;
    public:
        EulerFlux(int,int,int,int,
                  int,int,int);
        ~EulerFlux();
        
    public:
        void combineFluxLR();
    };
    
    class EulerData3D
    {
        friend class physicsModel;
    protected:
        Field jacobian;
        
        VectorField xi_xyz;
        VectorField eta_xyz;
        VectorField zeta_xyz;
        
        VectorField W_Euler;//primitive values (rho,u,v,w,p,...)
        
        VectorField W0_Euler;//acoustics values (T,e,alpha)
        
        //Intermediate fields for Reimann problem
        EulerFlux Flux_xi;
        EulerFlux Flux_eta;
        EulerFlux Flux_zeta;
        
        //invicid derivatives
        VectorField dfdxi;
        VectorField dgdeta;
        VectorField dhdzeta;
        
        double dt;
    public:
        EulerData3D( int nx, int ny, int nz, int addEq);
        ~EulerData3D();
        
        //used by field operators and Riemann solvers
        EulerFlux& getFluxXi();
        EulerFlux& getFluxEta();
        EulerFlux& getFluxZeta();
        
        //used by field operators and PDE
        virtual VectorField& getDrivativeXi();
        virtual VectorField& getDrivativeEta();
        virtual VectorField& getDrivativeZeta();
        
        void setDerivativesInv();
        
        void setDerivativesXiInv();
        void setDerivativesEtaInv();
        void setDerivativesZetaInv();
        
        //used by Riemann solvers
        VectorField& getPrimatives();
        VectorField& getAcoustics();
        
        virtual void solve(PDEData3d &,
                           fieldOperator3d &,
                           bufferData &,
                           physicsModel &);
    protected:
        void solveRiemann(PDEData3d &,
                          physicsModel &);
        
        void solveInv(fieldOperator3d &,
                      bufferData &);
        
        virtual void solveVis(fieldOperator3d &,
                              physicsModel &,
                              bufferData &);
    private:
        virtual void solveLocal();
    };
    
    class PDEData3d
    {
        int nEquations;
        
        //a PDE is consist of such four vectors
        // other vectors are intermediate data
        /*
         dQ   df     dg     dh
         -- + ---- + ---- + ----- = source
         dt   dxi    deta   dzeta
         */
        
        //current vector
        VectorField Q_Euler;
        
        VectorField dfdxi;
        VectorField dgdeta;
        VectorField dhdzeta;
        
        VectorField RHS;
        
    public:
        PDEData3d();
        PDEData3d(int,int,int,int);
        ~PDEData3d();
        
        
        void initPDEData3d(int,int,int,int);
        
        VectorField& getRHS(EulerData3D &);
        VectorField& getQ();
        virtual void solveInv(fieldOperator3d &op);
        
    protected:
        void setDrivatives(EulerData3D &);
        void setRHS(EulerData3D &);
    };
    
}

#endif