#ifndef Euler3d_h
#define Euler3d_h
#include <cstdlib>
#include "fieldOperator.h"
#include "MPICommunicator.h"

namespace nuc3d
{
    class EulerData3D;
    class EulerFlux;
    class PDEData3d;
    class physicsModel;
    
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
        
        //used by Riemann solvers
        VectorField& getPrimatives();
        VectorField& getAcoustics();
        
        
        
        virtual void solve(PDEData3d &,
                           fieldOperator3d &,
                           std::vector<bufferData> &,
                           physicsModel &,
                           MPIComunicator3d_nonblocking &);
    protected:
        void solveRiemann(PDEData3d &,
                          physicsModel &);
        
        void solveInv(fieldOperator3d &,
                      std::vector<bufferData> &,
                      MPIComunicator3d_nonblocking &);
    private:
        virtual void solveRHS(PDEData3d &);
        
        void solveInvicidFluxL(EulerFlux &,
                               fieldOperator3d &myOP,
                               std::vector<bufferData> &,
                               MPIComunicator3d_nonblocking &,
                               int );
        
        void solveInvicidFluxR(EulerFlux &,
                               fieldOperator3d &myOP,
                               std::vector<bufferData> &,
                               MPIComunicator3d_nonblocking &,
                               int );
        
        void setDerivativesInv();
        
        void setDerivativesXiInv();
        void setDerivativesEtaInv();
        void setDerivativesZetaInv();
    };
    }

#endif