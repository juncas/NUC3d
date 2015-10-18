#ifndef Euler3d_h
#define Euler3d_h

#include <cstdlib>
#include <iostream>
#include "field.h"


namespace nuc3d
{
    class PDEData3d
    {
        friend class MeshBlock;
        friend class meshGrid;
        int nEquations;
        
        //a PDE is consist of such four vectors
        // other vectors are intermediate data
        /*
         dQ   df     dg     dh
         -- + ---- + ---- + ----- = 0
         dt   dxi    deta   dzeta
         */
        
        //current vector
        VectorField Q_Euler;
        
        //working vector
        VectorField Q0_Euler;
        
        
        VectorField dfdxi;
        VectorField dgdeta;
        VectorField dhdzeta;
        
    public:
        PDEData3d(int,int,int,int);
        ~PDEData3d();
    };
    
    class EulerFlux
    {
        VectorField FluxL;
        VectorField FluxR;
        VectorField reconstFluxL;
        VectorField reconstFluxR;
        VectorField reconstFlux;
    public:
        EulerFlux(int,int,int,int);
        ~EulerFlux();
        
    public:
        void combineFluxLR();
    };
    
    class EulerData3D
    {
        friend class MeshBlock;
        friend class meshGrid;
        friend class physicsModel;
    protected:
        Field jacobian;
        
        VectorField xi_xyz;
        VectorField eta_xyz;
        VectorField zeta_xyz;
        
        
        VectorField W_Euler;//primitive values (rho,u,v,w,p)
        
        VectorField W0_Euler;//primitive values (T,e,alpha)
        
        
        //Intermediate fields for Reimann problem
        EulerFlux Flux_xi;
        EulerFlux Flux_eta;
        EulerFlux Flux_zeta;
        
        VectorField RHS;
        
        double dt;
        double maxEigen_xi;
        double maxEigen_eta;
        double maxEigen_zeta;
        
        
    public:
        EulerData3D( int nx, int ny, int nz, int addEq);
        ~EulerData3D();
    };
    
    class EulerReactiveData3D : virtual public EulerData3D
    {
        friend class MeshBlock;
        friend class meshGrid;
        friend class physicsModel;
        
        VectorField sourceTerm;
        
    public:
        EulerReactiveData3D();
        ~EulerReactiveData3D();
        
    };
    
    class NaiverStokesData3d : virtual public EulerData3D
    {
        friend class MeshBlock;
        friend class meshGrid;
        friend class physicsModel;
        
    protected:
        VectorField du;
        VectorField dv;
        VectorField dw;
        VectorField dt;
        
        //viscous stresses
        VectorField tau;
        
        //viscous fluxes
        VectorField fv;
        VectorField gv;
        VectorField hv;
        
        //viscous derivatives
        VectorField dfvdxi;
        VectorField dgvdeta;
        VectorField dhvdzeta;
        
    public:
        NaiverStokesData3d();
        ~NaiverStokesData3d();
        
    public:
        void getViscousStress();
        void getViscousFluxex();
        
    };
    
    class NaiverStokesReactiveData3d : public EulerReactiveData3D,public NaiverStokesData3d
    {
        
    public:
        NaiverStokesReactiveData3d();
        ~NaiverStokesReactiveData3d();
    };
}

#endif