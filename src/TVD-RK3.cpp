#ifndef tvdrk3_cpp
#define tvdrk3_cpp
#include "TVD-RK3.h"

nuc3d::tvdrk3rd::tvdrk3rd()
{
}


void nuc3d::tvdrk3rd::initial(const VectorField &U)
{
    for(int i=0;i<3;i++)
        ui.push_back(U);
    u=U;
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


void nuc3d::tvdrk3rd::integrationAll(const VectorField &RHS, // Right-hand-side: l*dt
                                           VectorField &Q, // u(nstep)
                                     double dt,
                                     int nstep)
{
    switch(nstep)
    {
        case 1:
            rk1st(RHS,Q,dt);
            break;
        case 2:
            rk2nd(RHS,Q,dt);
            break;
        case 3:
            rk3rd(RHS,Q,dt);
            break;
    }
    
}

void nuc3d::tvdrk3rd::rk1st(const VectorField &RHS,
                                  VectorField &un,
                                    double dt)
{
    int nx=u0.begin()->getSizeX();
    int ny=u0.begin()->getSizeY();
    int nz=u0.begin()->getSizeZ();
    
    
    for (auto iter=ui[0].begin(); iter!=ui[0].end(); ++iter)
    {
        u0=un;
        for(int k=0;k<nz;++k)
        {
            for(int j=0;j<ny;++j)
            {
                for(int i=0;i<nx;++i)
                {
                    double rhs=-RHS[iter-ui[0].begin()].getValue(i,j,k);
                    double q=iter->getValue(i, j, k);
                    
                    double u1 = coeff_tvdrk3_alpha0[nstep][0]*q
                    +coeff_tvdrk3_alpha0[nstep][1]*rhs*dt;
                    
                    iter->setValue(i,j,k,u1);
                }
            }
        }
        un=ui[0];
    }
}

void nuc3d::tvdrk3rd::rk2nd(const VectorField &RHS,
                            VectorField &un,
                            double dt)
{
    int nx=u0.begin()->getSizeX();
    int ny=u0.begin()->getSizeY();
    int nz=u0.begin()->getSizeZ();
    
    
    for (auto iter=ui[1].begin(); iter!=ui[1].end(); ++iter)
    {
        for(int k=0;k<nz;++k)
        {
            for(int j=0;j<ny;++j)
            {
                for(int i=0;i<nx;++i)
                {
                    double rhs=-RHS[iter-ui[1].begin()].getValue(i,j,k);
                    double q0=u0[iter-ui[1].begin()].getValue(i, j, k);
                    double q1=ui[0][iter-ui[1].begin()].getValue(i, j, k);
                    
                    double u1 = coeff_tvdrk3_alpha0[1][0]*q0
                    +coeff_tvdrk3_alpha0[1][1]*(rhs*dt+q1);
                    
                    iter->setValue(i,j,k,u1);
                }
            }
        }
        un=ui[1];
    }
}

void nuc3d::tvdrk3rd::rk3rd(const VectorField &RHS,
                            VectorField &un,
                            double dt)
{
    int nx=u0.begin()->getSizeX();
    int ny=u0.begin()->getSizeY();
    int nz=u0.begin()->getSizeZ();
    
    
    for (auto iter=ui[2].begin(); iter!=ui[2].end(); ++iter)
    {
        for(int k=0;k<nz;++k)
        {
            for(int j=0;j<ny;++j)
            {
                for(int i=0;i<nx;++i)
                {
                    double rhs=-RHS[iter-ui[2].begin()].getValue(i,j,k);
                    double q0=u0[iter-ui[2].begin()].getValue(i, j, k);
                    double q2=ui[1][iter-ui[2].begin()].getValue(i, j, k);

                    
                    double q1 = coeff_tvdrk3_alpha0[2][0]*q0
                    +coeff_tvdrk3_alpha0[2][1]*(rhs*dt+q2);
                    
                    iter->setValue(i,j,k,q1);
                }
            }
        }
        un=ui[0];
    }
}

#endif