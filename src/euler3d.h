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
        VectorField Q_Euler;
        
        VectorField dfdxi;
        VectorField dgdeta;
        VectorField dhdzeta;
    
        
        VectorField Q0_Euler;
        
        //temperary conservative fields for time marching
        VectorField Qi_Euler;
        VectorField Qii_Euler;
        VectorField Qiii_Euler;
    public:
        PDEData3d(int,int,int,int);
        ~PDEData3d();
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

		VectorField W0_Euler;//primitive values (T,e,alpha)

		VectorField W_Euler;//primitive values (rho,u,v,w,p)
        

		//Intermediate fields for Reimann problem
		VectorField FluxL_xi;
		VectorField FluxR_xi;
		VectorField reconstFluxL_xi;
		VectorField reconstFluxR_xi;
		VectorField reconstFlux_xi;

		VectorField FluxL_eta;
		VectorField FluxR_eta;
		VectorField reconstFluxL_eta;
		VectorField reconstFluxR_eta;
		VectorField reconstFlux_eta;

		VectorField FluxL_zeta;
		VectorField FluxR_zeta;
		VectorField reconstFluxL_zeta;
		VectorField reconstFluxR_zeta;
		VectorField reconstFlux_zeta;
		VectorField RHS;

		double dt;
		double maxEigen_xi;
		double maxEigen_eta;
		double maxEigen_zeta;


	public:
		EulerData3D( int nx, int ny, int nz, int addEq);
		~EulerData3D();
    };
    
    class EulerReactiveData3D : public EulerData3D
    {
        friend class MeshBlock;
        friend class meshGrid;
        friend class physicsModel;
        
        VectorField sourceTerm;
        
        
        
    };

    class NaiverStokesData3d : public EulerData3D
    {
        friend class MeshBlock;
        friend class meshGrid;
        friend class physicsModel;
        
    protected:
        VectorField du;
        VectorField dv;
        VectorField dw;
        VectorField dt;
        
        VectorField tau;
        
        VectorField fv;
        VectorField gv;
        VectorField hv;
        
        
    
        
        
        
        
    };
    
    class NaiverStokesReactiveData3d : public NaiverStokesData3d
    {
        VectorField sourceTerm;

        
    };
}

#endif