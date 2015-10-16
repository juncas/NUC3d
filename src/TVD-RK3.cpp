#ifndef tvdrk3_cpp
#define tvdrk3_cpp
#include "TVD-RK3.h"

nuc3d::tvdrk3::tvdrk3()
{

}


nuc3d::tvdrk3::~tvdrk3()
{
	
}

void nuc3d::tvdrk3::getRHS(
                const Field & fieldIN0, // dfdx
                const Field & fieldIN1, // dgdy
                const Field & fieldIN2, // dhdz
                const double dt,
                Field & fieldOUT)  // Right-hand-side
{
 	int nx=fieldOUT.getSizeX();
    int ny=fieldOUT.getSizeY();
    int nz=fieldOUT.getSizeZ();

    for(int k=0;k<nz;++k)
    {
        for(int j=0;j<ny;++j)
        { 
            for(int i=0;i<nx;++i)
            {
            	double dfdx=fieldIN0.getValue(i,j,k);
            	double dgdy=fieldIN1.getValue(i,j,k);
            	double dhdz=fieldIN2.getValue(i,j,k);

            	double rhs=-(dfdx+dgdy+dhdz)*dt;

           		fieldOUT.setValue(i,j,k,rhs);
            }
    	}
	}


}


void nuc3d::tvdrk3::integrationAll(
                const Field &RHS, // Right-hand-side: l*dt
                const Field &fieldIN0, // u(nstep)
                const Field &fieldIN1, // u_n
                int nstep, // n th step
                Field &fieldOUT) // output: u(nstep+1)
{
 	int nx=RHS.getSizeX();
    int ny=RHS.getSizeY();
    int nz=RHS.getSizeZ();


    if( (nx==(fieldOUT.getSizeX()))&&
        (ny==(fieldOUT.getSizeY()))&&
        (nz==(fieldOUT.getSizeZ())) )
    {
        for(int k=0;k<nz;++k)
        {
            for(int j=0;j<ny;++j)
            { 
                for(int i=0;i<nx;++i)
                {
                	double l=RHS.getValue(i,j,k);
                	double u_step=fieldIN0.getValue(i,j,k);
                	double u_n0=fieldIN1.getValue(i,j,k);


                	double u_n = coeff_tvdrk3_alpha0[nstep][0]*u_n0
                			    +coeff_tvdrk3_alpha0[nstep][1]*(l+u_step);

               		fieldOUT.setValue(i,j,k,u_n);
                }
        	}
    	}

    }

}

#endif