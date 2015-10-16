#ifndef Euler3d_h
#define Euler3d_h

#include <cstdlib>
#include <iostream>
#include "field.h" 


namespace nuc3d 
{

	class EulerData3D
	{
		friend class MeshBlock;
		friend class meshGrid;
		friend class physicsModel;
		
		int nEquations;
		
		Field jacobian;

		VectorField xi_xyz;
		VectorField eta_xyz;
		VectorField zeta_xyz;

		VectorField W0_Euler;//primitive values (T,e,alpha)

	//flow fields;

		VectorField W_Euler;//primitive values (rho,u,v,w,p)

		//conservative values
		VectorField Q_Euler;

		//temperary conservative fields for time marching
		VectorField Qi_Euler;
		VectorField Qii_Euler;
		VectorField Qiii_Euler;

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

		VectorField dfdxi;
		VectorField dfdeta;
		VectorField dfdzeta;
		VectorField RHS;

		double dt;
		double maxEigen_xi;
		double maxEigen_eta;
		double maxEigen_zeta;


	public:
		EulerData3D(const int nx,const int ny,const int nz,const int addEq);
		~EulerData3D();
	public:


	private:
		
	};
}

#endif