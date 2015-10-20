#ifndef Euler3d_h
#define Euler3d_h

#include <cstdlib>
#include <iostream>
#include "physicsModel.h"

namespace nuc3d
{
    class meshGrid;
    class MeshBlock;
    
    class EulerFlux
    {
        friend class physicsModel;
        
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
        
        virtual void solveLocal();
    };
    
    class EulerReactiveData3D : virtual public EulerData3D
    {
        friend class physicsModel;
        
    protected:
        VectorField source_xi;
        VectorField source_eta;
        VectorField source_zeta;
        
    public:
        EulerReactiveData3D(int,int,int,int);
        ~EulerReactiveData3D();
        
        virtual VectorField& getDrivativeXi();
        virtual VectorField& getDrivativeEta();
        virtual VectorField& getDrivativeZeta();
        virtual void solveLocal();
    protected:
        void setSource();
    };
    
    class NaiverStokesData3d : virtual public EulerData3D
    {
        friend class physicsModel;
        
    protected:
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
        
    protected:
        void setViscousFluxes();
    };
    
    class NaiverStokesReactiveData3d : public EulerReactiveData3D, public NaiverStokesData3d
    {
    public:
        NaiverStokesReactiveData3d(int,int,int,int);
        ~NaiverStokesReactiveData3d();
        
        virtual VectorField& getDrivativeXi();
        virtual VectorField& getDrivativeEta();
        virtual VectorField& getDrivativeZeta();
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
        
    protected:
        void setDrivatives(EulerData3D &);
        void setRHS(EulerData3D &);
    };

}

#endif