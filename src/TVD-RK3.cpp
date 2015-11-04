#ifndef tvdrk3_cpp
#define tvdrk3_cpp
#include "TVD-RK3.h"

nuc3d::tvdrk3rd::tvdrk3rd()
{
}


nuc3d::tvdrk3rd::~tvdrk3rd()
{
    
}

void nuc3d::tvdrk3rd::integrationAll(const VectorField &RHS, // Right-hand-side: l*dt
                                     const VectorField &Q, // u(nstep)
                                     VectorField &Qi, //u_i
                                     double dt,
                                     int nstep)
{
    switch(nstep)
    {
        case 1:
            rk1st(RHS,Q,Qi,dt);
            break;
        case 2:
            rk2nd(RHS,Q,Qi,dt);
            break;
        case 3:
            rk3rd(RHS,Q,Qi,dt);
            break;
    }
    
}

void nuc3d::tvdrk3rd::rk1st(const VectorField &RHS,
                            const VectorField &u_n, // u_n
                            VectorField &u_i,
                            double dt)
{
    int nx=u_n.begin()->getSizeX();
    int ny=u_n.begin()->getSizeY();
    int nz=u_n.begin()->getSizeZ();
    
    auto beg=u_i.begin();
    auto end=u_i.end();
    
    for (auto iter=beg; iter!=end; ++iter)
    {
        for(int k=0;k<nz;++k)
        {
            for(int j=0;j<ny;++j)
            {
                for(int i=0;i<nx;++i)
                {
                    double rhs=-RHS[iter-beg].getValue(i,j,k);
                    double u=iter->getValue(i, j, k);
                    
                    double u1 = coeff_tvdrk3_alpha0[0][0]*u
                    +coeff_tvdrk3_alpha0[0][1]*rhs*dt;
                    
                    iter->setValue(i,j,k,u1);
                }
            }
        }
    }
}

void nuc3d::tvdrk3rd::rk2nd(const VectorField &RHS,
                            const VectorField &u_n, // u_n
                            VectorField &u_i,
                            double dt)
{
    int nx=u_n.begin()->getSizeX();
    int ny=u_n.begin()->getSizeY();
    int nz=u_n.begin()->getSizeZ();
    
    auto beg=u_i.begin();
    auto end=u_i.end();
    
    
    for (auto iter=beg; iter!=end; ++iter)
    {
        for(int k=0;k<nz;++k)
        {
            for(int j=0;j<ny;++j)
            {
                for(int i=0;i<nx;++i)
                {
                    double rhs=-RHS[iter-beg].getValue(i,j,k);
                    double u=u_n[iter-beg].getValue(i, j, k);
                    double u1=iter->getValue(i, j, k);
                    
                    double u2 = coeff_tvdrk3_alpha0[1][0]*u
                    +coeff_tvdrk3_alpha0[1][1]*(rhs*dt+u1);
                    
                    iter->setValue(i,j,k,u2);
                }
            }
        }
    }
}

void nuc3d::tvdrk3rd::rk3rd(const VectorField &RHS,
                            const VectorField &u_n, // u_n
                            VectorField &u_i,
                            double dt)
{
    int nx=u_n.begin()->getSizeX();
    int ny=u_n.begin()->getSizeY();
    int nz=u_n.begin()->getSizeZ();
    
    auto beg=u_i.begin();
    auto end=u_i.end();
    
    
    for (auto iter=beg; iter!=end; ++iter)
    {
        for(int k=0;k<nz;++k)
        {
            for(int j=0;j<ny;++j)
            {
                for(int i=0;i<nx;++i)
                {
                    double rhs=-RHS[iter-beg].getValue(i,j,k);
                    double u=u_n[iter-beg].getValue(i, j, k);
                    double u2=iter->getValue(i, j, k);
                    
                    double u3 = coeff_tvdrk3_alpha0[2][0]*u
                    +coeff_tvdrk3_alpha0[2][1]*(rhs*dt+u2);
                    
                    iter->setValue(i,j,k,u3);
                }
            }
        }
    }
}

#endif