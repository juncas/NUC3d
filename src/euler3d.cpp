#ifndef euler3d_cpp
#define euler3d_cpp
#include "euler3d.h"
/**************************************************************************************
 Member functions of class: PDEData3d
 **************************************************************************************/
nuc3d::PDEData3d::PDEData3d(int nx0,int ny0,int nz0,int neqs):
nEquations(neqs),
Q_Euler(neqs,Field(nx0,ny0,nz0)),
dfdxi(neqs,Field(nx0,ny0,nz0)),
dgdeta(neqs,Field(nx0,ny0,nz0)),
dhdzeta(neqs,Field(nx0,ny0,nz0))
{}


void nuc3d::PDEData3d::setRHS(EulerData3D &myFluxes)
{
    setDrivatives(myFluxes);
    for(auto iter=RHS.begin();iter!=RHS.end();iter++)
    {
        *iter=dfdxi[iter-RHS.begin()];
        *iter+=dgdeta[iter-RHS.begin()];
        *iter+=dhdzeta[iter-RHS.begin()];
    }
}

void nuc3d::PDEData3d::setDrivatives(EulerData3D &myFluxes)
{
    dfdxi=myFluxes.getDrivativeXi();
    dfdxi=myFluxes.getDrivativeEta();
    dfdxi=myFluxes.getDrivativeZeta();
}

nuc3d::VectorField& nuc3d::PDEData3d::getRHS(EulerData3D &myFluxes)
{
    setRHS(myFluxes);
    return RHS;
}

nuc3d::PDEData3d::~PDEData3d()
{}
/**************************************************************************************
 Member functions of class: EulerFlux
 **************************************************************************************/
nuc3d::EulerFlux::EulerFlux(int nx0,int ny0,int nz0,int neqs,
                            int xdir,int ydir,int zdir):
FluxL(neqs,Field(nx0,ny0,nz0)),
FluxR(neqs,Field(nx0,ny0,nz0)),
reconstFluxL(neqs,Field(nx0+xdir,ny0+ydir,nz0+zdir)),
reconstFluxR(neqs,Field(nx0+xdir,ny0+ydir,nz0+zdir)),
reconstFlux(neqs,Field(nx0,ny0,nz0))
{
    
}

void nuc3d::EulerFlux::combineFluxLR()
{
    for(auto iter=reconstFlux.begin();iter!=reconstFlux.end();iter++)
        *iter=reconstFluxL[iter-reconstFlux.begin()]
        +reconstFluxR[iter-reconstFlux.begin()];
}

nuc3d::EulerFlux::~EulerFlux()
{}

/**************************************************************************************
 Member functions of class: EulerData3D
 **************************************************************************************/
nuc3d::EulerData3D::EulerData3D( int nx0, int ny0, int nz0, int neqs):
jacobian(Field(nx0,ny0,nz0)),
xi_xyz(3,Field(nx0,ny0,nz0)),
eta_xyz(3,Field(nx0,ny0,nz0)),
zeta_xyz(3,Field(nx0,ny0,nz0)),
W0_Euler(3,Field(nx0,ny0,nz0)),
W_Euler(neqs,Field(nx0,ny0,nz0)),
Flux_xi(nx0,ny0,nz0,neqs,1,0,0),
Flux_eta(nx0,ny0,nz0,neqs,0,1,0),
Flux_zeta(nx0,ny0,nz0,neqs,0,0,1),
dfdxi(neqs,Field(nx0,ny0,nz0)),
dgdeta(neqs,Field(nx0,ny0,nz0)),
dhdzeta(neqs,Field(nx0,ny0,nz0)),
dt(0.0),
maxEigen_xi(0.0),
maxEigen_eta(0.0),
maxEigen_zeta(0.0)
{
}


nuc3d::VectorField& nuc3d::EulerData3D::getDrivativeXi()
{
    return this->dfdxi;
}

nuc3d::VectorField& nuc3d::EulerData3D::getDrivativeEta()
{
    return this->dgdeta;
}

nuc3d::VectorField& nuc3d::EulerData3D::getDrivativeZeta()
{
    return this->dhdzeta;
}

nuc3d::EulerFlux& nuc3d::EulerData3D::getFluxXi()
{
    return Flux_xi;
}

nuc3d::EulerFlux& nuc3d::EulerData3D::getFluxEta()
{
    return Flux_eta;
}

nuc3d::EulerFlux& nuc3d::EulerData3D::getFluxZeta()
{
    return Flux_zeta;
}

nuc3d::VectorField& nuc3d::EulerData3D::getPrimatives()
{
    return W_Euler;
}


nuc3d::VectorField& nuc3d::EulerData3D::getAcoustics()
{
    return W0_Euler;
}


nuc3d::EulerData3D::~EulerData3D()
{}
/**************************************************************************************
 Member functions of class: EulerReactiveData3D
 **************************************************************************************/
nuc3d::EulerReactiveData3D::EulerReactiveData3D( int nx0, int ny0, int nz0, int neqs):
EulerData3D(nx0,ny0,nz0,neqs),
source_xi(neqs,Field(nx0,ny0,nz0,0.0)),
source_eta(neqs,Field(nx0,ny0,nz0,0.0)),
source_zeta(neqs,Field(nx0,ny0,nz0,0.0))
{
    
}


nuc3d::VectorField& nuc3d::EulerReactiveData3D::getDrivativeXi()
{
    for(auto iter=dfdxi.begin();iter!=dfdxi.end();iter++)
        *iter-=source_xi[iter-dfdxi.begin()];
    
    return this->dfdxi;
}

nuc3d::VectorField& nuc3d::EulerReactiveData3D::getDrivativeEta()
{
    for(auto iter=dgdeta.begin();iter!=dgdeta.end();iter++)
        *iter-=source_eta[iter-dfdxi.begin()];
    
    return this->dgdeta;
}

nuc3d::VectorField& nuc3d::EulerReactiveData3D::getDrivativeZeta()
{
    for(auto iter=dhdzeta.begin();iter!=dhdzeta.end();iter++)
        *iter-=source_zeta[iter-dfdxi.begin()];
    
    return this->dhdzeta;
}

void nuc3d::EulerReactiveData3D::solveLocal()
{
    setSource();
}

void nuc3d::EulerReactiveData3D::setSource()
{
    
}

nuc3d::EulerReactiveData3D::~EulerReactiveData3D()
{
    
}

/**************************************************************************************
 Member functions of class: NaiverStokesData3d
 **************************************************************************************/
nuc3d::NaiverStokesData3d::NaiverStokesData3d(int nx0, int ny0, int nz0, int neqs):
miu(nx0,ny0,nz0,1.0),
coeff(nx0,ny0,nz0,1.0),
EulerData3D(nx0,ny0,nz0,neqs),
du(3,Field(nx0,ny0,nz0)),
dv(3,Field(nx0,ny0,nz0)),
dw(3,Field(nx0,ny0,nz0)),
dT(3,Field(nx0,ny0,nz0)),
tau(9,Field(nx0,ny0,nz0)),
Flux_xi_vis(neqs,Field(nx0,ny0,nz0)),
Flux_eta_vis(neqs,Field(nx0,ny0,nz0)),
Flux_zeta_vis(neqs,Field(nx0,ny0,nz0)),
dfvdxi(neqs,Field(nx0,ny0,nz0)),
dgvdeta(neqs,Field(nx0,ny0,nz0)),
dhvdzeta(neqs,Field(nx0,ny0,nz0))
{
    
}

void nuc3d::NaiverStokesData3d::solveLocal()
{
    setViscousFluxes();
}

void nuc3d::NaiverStokesData3d::setViscousFluxes()
{
    int nx=miu.getSizeX();
    int ny=miu.getSizeY();
    int nz=miu.getSizeZ();
    
    for (int k=0; k<nz; k++) {
        for (int j=0; j<ny; j++) {
            for (int i=0; i<nx; i++) {
                double u=W_Euler[1].getValue(i, j, k);
                double v=W_Euler[2].getValue(i, j, k);
                double w=W_Euler[3].getValue(i, j, k);
                
                double uxi=du[0].getValue(i, j, k);
                double ueta=du[1].getValue(i, j, k);
                double uzeta=du[2].getValue(i, j, k);
                
                double vxi=dv[0].getValue(i, j, k);
                double veta=dv[1].getValue(i, j, k);
                double vzeta=dv[2].getValue(i, j, k);
                
                double wxi=dw[0].getValue(i, j, k);
                double weta=dw[1].getValue(i, j, k);
                double wzeta=dw[2].getValue(i, j, k);
                
                double txi=dT[0].getValue(i, j, k);
                double teta=dT[1].getValue(i, j, k);
                double tzeta=dT[2].getValue(i, j, k);
                
                double xi_x=xi_xyz[0].getValue(i, j, k);
                double eta_x=xi_xyz[1].getValue(i, j, k);
                double zeta_x=xi_xyz[2].getValue(i, j, k);
                
                double xi_y=eta_xyz[0].getValue(i, j, k);
                double eta_y=eta_xyz[1].getValue(i, j, k);
                double zeta_y=eta_xyz[2].getValue(i, j, k);
                
                double xi_z=zeta_xyz[0].getValue(i, j, k);
                double eta_z=zeta_xyz[1].getValue(i, j, k);
                double zeta_z=zeta_xyz[2].getValue(i, j, k);
                
                double jac=jacobian.getValue(i, j, k);
                
                double miu0=miu.getValue(i, j, k);
                double coeff0=coeff.getValue(i, j, k);
                
                double ux=uxi*xi_x+ueta*eta_x+uzeta*zeta_x;
                double uy=uxi*xi_y+ueta*eta_y+uzeta*zeta_y;
                double uz=uxi*xi_z+ueta*eta_z+uzeta*zeta_z;
                
                double vx=vxi*xi_x+veta*eta_x+vzeta*zeta_x;
                double vy=vxi*xi_y+veta*eta_y+vzeta*zeta_y;
                double vz=vxi*xi_z+veta*eta_z+vzeta*zeta_z;
                
                double wx=wxi*xi_x+weta*eta_x+wzeta*zeta_x;
                double wy=wxi*xi_y+weta*eta_y+wzeta*zeta_y;
                double wz=wxi*xi_z+weta*eta_z+wzeta*zeta_z;
                
                double grad=ux+vy+wz;
                
                double tau_xx=miu0*(2.0*ux-2.0/3.0*grad);
                double tau_xy=miu0*(uy+vx-2.0/3.0*grad);
                double tau_xz=miu0*(uz+wx-2.0/3.0*grad);
                double tau_yy=miu0*(2.0*vy-2.0/3.0*grad);
                double tau_yz=miu0*(vz+wy-2.0/3.0*grad);
                double tau_zz=miu0*(2.0*wz-2.0/3.0*grad);
                
                double tau_tx=coeff0*(txi*xi_x+teta*eta_x+tzeta*zeta_x);
                double tau_ty=coeff0*(txi*xi_y+teta*eta_y+tzeta*zeta_y);
                double tau_tz=coeff0*(txi*xi_z+teta*eta_z+tzeta*zeta_z);
                
                tau[0].setValue(i,j, k, tau_xx);
                tau[1].setValue(i,j, k, tau_xy);
                tau[2].setValue(i,j, k, tau_xz);
                tau[3].setValue(i,j, k, tau_yy);
                tau[4].setValue(i,j, k, tau_yz);
                tau[5].setValue(i,j, k, tau_zz);
                tau[6].setValue(i,j, k, tau_tx);
                tau[7].setValue(i,j, k, tau_ty);
                tau[8].setValue(i,j, k, tau_tz);
                
                double fv[5];
                double gv[5];
                double hv[5];
                
                fv[0]=0.0;
                fv[1]=(tau_xx*xi_x+tau_xy*xi_y+tau_xz*xi_z)/jac;
                fv[2]=(tau_xz*xi_x+tau_yz*xi_y+tau_zz*xi_z)/jac;
                fv[3]=(tau_xx*xi_x+tau_xy*xi_y+tau_xz*xi_z)/jac;
                fv[4]=(u*tau_xx-tau_tx)*xi_x
                +(v*tau_xy-tau_ty)*xi_y
                +(w*tau_xz-tau_tz)*xi_z;
                
                gv[0]=0.0;
                gv[1]=(tau_xx*eta_x+tau_xy*eta_y+tau_xz*eta_z)/jac;
                gv[2]=(tau_xz*eta_x+tau_yz*eta_y+tau_zz*eta_z)/jac;
                gv[3]=(tau_xx*eta_x+tau_xy*eta_y+tau_xz*eta_z)/jac;
                gv[4]=(u*tau_xx-tau_tx)*eta_x
                +(v*tau_xy-tau_ty)*eta_y
                +(w*tau_xz-tau_tz)*eta_z;

                hv[0]=0.0;
                hv[1]=(tau_xx*zeta_x+tau_xy*zeta_y+tau_xz*zeta_z)/jac;
                hv[2]=(tau_xz*zeta_x+tau_yz*zeta_y+tau_zz*zeta_z)/jac;
                hv[3]=(tau_xx*zeta_x+tau_xy*zeta_y+tau_xz*zeta_z)/jac;
                hv[4]=(u*tau_xx-tau_tx)*zeta_x
                +(v*tau_xy-tau_ty)*zeta_y
                +(w*tau_xz-tau_tz)*zeta_z;
                
                Flux_xi_vis[0].setValue(i, j, k,fv[0]);
                Flux_xi_vis[1].setValue(i, j, k,fv[1]);
                Flux_xi_vis[2].setValue(i, j, k,fv[2]);
                Flux_xi_vis[3].setValue(i, j, k,fv[3]);
                Flux_xi_vis[4].setValue(i, j, k,fv[4]);
                for (auto iter=Flux_xi_vis.begin()+5; iter!=Flux_xi_vis.end(); ++iter) {
                    iter->setValue(i, j, k, 0.0);
                }
                
                Flux_eta_vis[0].setValue(i, j, k,fv[0]);
                Flux_eta_vis[1].setValue(i, j, k,fv[1]);
                Flux_eta_vis[2].setValue(i, j, k,fv[2]);
                Flux_eta_vis[3].setValue(i, j, k,fv[3]);
                Flux_eta_vis[4].setValue(i, j, k,fv[4]);
                for (auto iter=Flux_eta_vis.begin()+5; iter!=Flux_eta_vis.end(); ++iter) {
                    iter->setValue(i, j, k, 0.0);
                }
                
                Flux_zeta_vis[0].setValue(i, j, k,fv[0]);
                Flux_zeta_vis[1].setValue(i, j, k,fv[1]);
                Flux_zeta_vis[2].setValue(i, j, k,fv[2]);
                Flux_zeta_vis[3].setValue(i, j, k,fv[3]);
                Flux_zeta_vis[4].setValue(i, j, k,fv[4]);
                for (auto iter=Flux_zeta_vis.begin()+5; iter!=Flux_zeta_vis.end(); ++iter) {
                    iter->setValue(i, j, k, 0.0);
                }
                
            }
        }
    }
    
}

nuc3d::VectorField& nuc3d::NaiverStokesData3d::getDrivativeXi()
{
    for(auto iter=dfdxi.begin();iter!=dfdxi.end();iter++)
        *iter-=dfvdxi[iter-dfdxi.begin()];
    
    return this->dfdxi;
}

nuc3d::VectorField& nuc3d::NaiverStokesData3d::getDrivativeEta()
{
    for(auto iter=dgdeta.begin();iter!=dgdeta.end();iter++)
        *iter-=dgdeta[iter-dfdxi.begin()];
    
    return this->dgdeta;
}

nuc3d::VectorField& nuc3d::NaiverStokesData3d::getDrivativeZeta()
{
    for(auto iter=dhdzeta.begin();iter!=dhdzeta.end();iter++)
        *iter-=dhvdzeta[iter-dfdxi.begin()];
    
    return this->dhdzeta;
}


nuc3d::NaiverStokesData3d::~NaiverStokesData3d()
{}
/**************************************************************************************
 Member functions of class: NaiverStokesReactiveData3d
 **************************************************************************************/
nuc3d::NaiverStokesReactiveData3d::NaiverStokesReactiveData3d(int nx0, int ny0, int nz0, int neqs):
EulerData3D(nx0,ny0,nz0,neqs),
EulerReactiveData3D(nx0,ny0,nz0,neqs),
NaiverStokesData3d(nx0,ny0,nz0,neqs)
{
    
}

nuc3d::VectorField& nuc3d::NaiverStokesReactiveData3d::getDrivativeXi()
{
    for(auto iter=dfdxi.begin();iter!=dfdxi.end();iter++)
    {
        *iter-=dfvdxi[iter-dfdxi.begin()];
        *iter-source_xi[iter-dfdxi.begin()];
    }
    
    return this->dfdxi;
}

nuc3d::VectorField& nuc3d::NaiverStokesReactiveData3d::getDrivativeEta()
{
    for(auto iter=dgdeta.begin();iter!=dgdeta.end();iter++)
    {
        *iter-=dgdeta[iter-dfdxi.begin()];
        *iter-=source_eta[iter-dfdxi.begin()];
    }
    
    return this->dgdeta;
}

nuc3d::VectorField& nuc3d::NaiverStokesReactiveData3d::getDrivativeZeta()
{
    for(auto iter=dhdzeta.begin();iter!=dhdzeta.end();iter++)
    {
        *iter-=dhvdzeta[iter-dfdxi.begin()];
        *iter-=source_zeta[iter-dfdxi.begin()];
    }
    
    return this->dhdzeta;
}

void nuc3d::NaiverStokesReactiveData3d::solveLocal()
{
    EulerReactiveData3D::solveLocal();
    NaiverStokesReactiveData3d::solveLocal();
}

nuc3d::NaiverStokesReactiveData3d::~NaiverStokesReactiveData3d()
{}
#endif