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
    class boundaryCondition;
    
    class EulerFlux
    {
        friend class physicsModel;
        friend class EulerData3D;
    public:
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
        int nx;
        int ny;
        int nz;
        EulerData3D( int, int , int , int );
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
        
        void initialxyz(VectorField &);
        Field& getJac(){return jacobian;};
        VectorField& getXi_xyz(){return xi_xyz;};
        VectorField& getEta_xyz(){return eta_xyz;};
        VectorField& getZeta_xyz(){return zeta_xyz;};
        
        virtual void solve(PDEData3d &,
                           fieldOperator3d &,
                           VectorBuffer &,
                           physicsModel &,
                           MPIComunicator3d_nonblocking &,
                           boundaryCondition &);
    protected:
        void solveRiemann(PDEData3d &,
                          physicsModel &);
        
        void solveCon2Prim(PDEData3d &myPDE,
                           physicsModel &myModel);
        
        void solveInv(fieldOperator3d &,
                      VectorBuffer &,
                      MPIComunicator3d_nonblocking &,
                      boundaryCondition &myBC);
        
        void setBoundaryCondition(PDEData3d &,
                                  physicsModel &myModel,
                                  std::vector<bufferData> &,
                                  boundaryCondition &);
        
        void getDt();
        
    private:
        void solveRHS(PDEData3d &);
        
        void solveInvicidFluxL(EulerFlux &,
                               fieldOperator3d &myOP,
                               VectorBuffer &,
                               MPIComunicator3d_nonblocking &,
                               boundaryCondition &myBC,
                               int );
        
        void solveInvicidFluxR(EulerFlux &,
                               fieldOperator3d &myOP,
                               VectorBuffer &,
                               MPIComunicator3d_nonblocking &,
                               boundaryCondition &myBC,
                               int );
        
        void setDerivativesInv();
        
        void setDerivativesXiInv();
        void setDerivativesEtaInv();
        void setDerivativesZetaInv();
    };
}

#endif