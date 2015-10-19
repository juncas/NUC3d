#ifndef tvdrk3_cpp
#define tvdrk3_cpp
#include "TVD-RK3.h"

nuc3d::tvdrk3rd::tvdrk3rd()
{
    
}


nuc3d::tvdrk3rd::~tvdrk3rd()
{
    
}

void nuc3d::tvdrk3rd::getRHS(
                             const VectorField & fieldIN0, // dfdx
                             const VectorField & fieldIN1, // dgdy
                             const VectorField & fieldIN2, // dhdz
                             const double dt,
                             VectorField & fieldOUT)  // Right-hand-side
{
    int nx=fieldOUT.begin()->getSizeX();
    int ny=fieldOUT.begin()->getSizeY();
    int nz=fieldOUT.begin()->getSizeZ();
    
    for(auto iter=fieldOUT.begin();iter!=fieldOUT.end();iter++)
    {
        for(int k=0;k<nz;++k)
        {
            for(int j=0;j<ny;++j)
            {
                for(int i=0;i<nx;++i)
                {
                    double dfdx=fieldIN0[iter-fieldOUT.begin()].getValue(i,j,k);
                    double dgdy=fieldIN1[iter-fieldOUT.begin()].getValue(i,j,k);
                    double dhdz=fieldIN2[iter-fieldOUT.begin()].getValue(i,j,k);
                    
                    double rhs=-(dfdx+dgdy+dhdz)*dt;
                    
                    iter->setValue(i,j,k,rhs);
                }
            }
        }
    }
    
}


void nuc3d::tvdrk3rd::integrationAll(
                                     const VectorField &RHS, // Right-hand-side: l*dt
                                     const VectorField &Ui, // u(nstep)
                                     const VectorField &U0, // u_n
                                     int nstep, // n th step
                                     VectorField &Uii) // output: u(nstep+1)
{
    int nx=RHS.begin()->getSizeX();
    int ny=RHS.begin()->getSizeY();
    int nz=RHS.begin()->getSizeZ();
    
    
    if( (nx==(Uii.begin()->getSizeX()))&&
        (ny==(Uii.begin()->getSizeY()))&&
        (nz==(Uii.begin()->getSizeZ())) )
    {
        for (auto iter=Uii.begin(); iter!=Uii.end(); ++iter)
        {
            for(int k=0;k<nz;++k)
            {
                for(int j=0;j<ny;++j)
                {
                    for(int i=0;i<nx;++i)
                    {
                        double l=-RHS[iter-Uii.begin()].getValue(i,j,k);
                        double u_step=Ui[iter-Uii.begin()].getValue(i,j,k);
                        double u_n0=U0[iter-Uii.begin()].getValue(i,j,k);
                        
                        
                        double u_n = coeff_tvdrk3_alpha0[nstep][0]*u_n0
                        +coeff_tvdrk3_alpha0[nstep][1]*(l+u_step);
                        
                        iter->setValue(i,j,k,u_n);
                    }
                }
            }

        }
        
    }
    
}

#endif