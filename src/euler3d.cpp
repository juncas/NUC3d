#ifndef euler3d_cpp
#define euler3d_cpp
#include "euler3d.h"
/**************************************************************************************
			Definition of class communicator: base class for communicators
**************************************************************************************/

nuc3d::EulerData3D::EulerData3D(const int nx,const int ny,const int nz,const int neqs):
dt(0.0),
maxEigen_xi(0.0),
maxEigen_eta(0.0),
maxEigen_zeta(0.0),
nEquations(neqs),
jacobian(Field(nx,ny,nz)),
xi_xyz(3,field3d<double>(nx,ny,nz)),
eta_xyz(3,field3d<double>(nx,ny,nz)),
zeta_xyz(3,field3d<double>(nx,ny,nz)),
W0_Euler(3,field3d<double>(nx,ny,nz)),
W_Euler(neqs,field3d<double>(nx,ny,nz)),
Q_Euler(neqs,field3d<double>(nx,ny,nz)),
Qi_Euler(neqs,field3d<double>(nx,ny,nz)),
Qii_Euler(neqs,field3d<double>(nx,ny,nz)),
Qiii_Euler(neqs,field3d<double>(nx,ny,nz)),
FluxL(neqs,field3d<double>(nx,ny,nz)),
FluxR(neqs,field3d<double>(nx,ny,nz)),
reconstFluxL(neqs,field3d<double>(nx,ny,nz)),
reconstFluxR(neqs,field3d<double>(nx,ny,nz)),
reconstFlux(neqs,field3d<double>(nx,ny,nz)),
dfdx(neqs,field3d<double>(nx,ny,nz)),
dfdy(neqs,field3d<double>(nx,ny,nz)),
dfdz(neqs,field3d<double>(nx,ny,nz)),
RHS(neqs,field3d<double>(nx,ny,nz)),
{
}

nuc3d::Euler3d::~Euler3d()
{

};
/**************************************************************************************
						Definition of constructors and destructors
**************************************************************************************/

/**************************************************************************************
						Definition of member functions
**************************************************************************************/
/**************************************************************************************
								End of definition
**************************************************************************************/


#endif