#ifndef physicsModel_cpp
#define physicsModel_cpp
#include "physicsModel.h"
#include "euler3d.h"
#include "PDEData3d.hpp"
#include <memory>
/**************************************************************************************
 Definition of class communicator: base class for communicators
 **************************************************************************************/
/**************************************************************************************
 Definition of constructors and destructors
 **************************************************************************************/
nuc3d::physicsModel::physicsModel() :
neqs(5),
myEosFWDMap(
            {
                {"IdealGas",&nuc3d::physicsModel::EoSIdealGasFWD},
                {"JWL",&nuc3d::physicsModel::EoSJWLFWD}
            }
            ),
myEosBWDMap(
            {
                {"IdealGas",&nuc3d::physicsModel::EoSIdealGasBWD},
                {"JWL",&nuc3d::physicsModel::EoSJWLBWD}
            }
            ),
myRiemannMap(
             {
                 { "AUSMp",&nuc3d::physicsModel::RiemannAUSMp },
                 { "AUSM",&nuc3d::physicsModel::RiemannAUSM },
                 { "LF",&nuc3d::physicsModel::RiemannLF}
             }
             ),
myRiemannPointMap(
                  {
                      { "AUSMp",&nuc3d::physicsModel::solveRiemannPointAUSMp},
                      { "AUSM",&nuc3d::physicsModel::solveRiemannPointAUSM},
                      { "LF",&::nuc3d::physicsModel::solveRiemannPointLF}
                  }
                  ),
myModelParameters(
                  {
                      {"Reynolds",1000},
                      {"Mach",0.5},
                      {"Pt",0.72},
                      {"Gamma",1.4},
                      {"T_ref",110.4},
                      {"T_inf",288.0},
                      {"T_wall",1.0}
                  }
                  ),
myVisModelMap(
{
    {"Sutherland",&nuc3d::physicsModel::sutherland},
    {"Constant",&nuc3d::physicsModel::constant}
}
              )
{
    std::string str;
    std::ifstream file("inp/PhysModel.in");
    //std::cout<<"Start reading PhysModel.in ..."<<std::endl;
    if (file)
    {
        // this order should not change
        file >> str >> neqs;
        file >> str >> myEoSName;
        file >> str >> myRiemannName;
        file >> str >> myModelName;
        file >> str >> myVisModelName;
        
        if(neqs<5)
        {
            std::cout<<"Equation number less than 5,using default"
            << std::endl;
            neqs=5;
        }
        
        
        if(myEosFWDMap.find(myEoSName)==myEosFWDMap.end())
        {
            std::cout<<"Equation of State name "<<myEoSName<<" does not exist ,using default!"
            << std::endl;
            myEoSName="IdealGas";
        }
        
        if(myRiemannMap.find(myRiemannName)==myRiemannMap.end())
        {
            std::cout<<"Riemann Solver name "<<myEoSName<<" does not exist ,using default!"
            << std::endl;
            myEoSName="AUSM";
        }
        
        
        if(myVisModelMap.find(myVisModelName)==myVisModelMap.end())
        {
            std::cout<<"Model name "<<myVisModelName<<" does not exist ,using default!"
            << std::endl;
            myVisModelName="Sutherland";
            
        }
        
        while (readPhysFile(file));
        //std::cout<<"PhysModel.in has been read!"<<std::endl;
    }
    else
    {
        std::cout << "IO file \'PhysModel.in\' does not exist!"
        << std::endl;
        exit(0);
    }
    file.close();
}

nuc3d::physicsModel::~physicsModel()
{
    
}
/**************************************************************************************
 Definition of member functions
 **************************************************************************************/
std::istream& nuc3d::physicsModel::readPhysFile(std::istream& ios)
{
    std::string word0;
    std::string word1;
    double value;
    ios >> word0 >> word1;
    std::istringstream value0(word1);
    value0 >> value;
    if (myModelParameters.find(word0) != myModelParameters.end())
    {
        myModelParameters[word0] = value;
    }
    else if(word0.size())
    {
        std::cout << "word " << word0 << " does not exist"<<std::endl;
        exit(0);
    }
    
    return ios;
    
}

void nuc3d::physicsModel::solve(PDEData3d &myPDE,EulerData3D *myEuler)
{
    con2prim(myEoSName,
             myEuler->jacobian,
             myPDE.getQ(),
             myEuler->W_Euler,
             myEuler->W0_Euler);
    
}

void nuc3d::physicsModel::getPrim(Field &jac,VectorField &Q,VectorField &Prim, VectorField &Acoust)
{
    con2prim(myEoSName,
             jac,
             Q,
             Prim,
             Acoust);
    
}

void  nuc3d::physicsModel::solveRiemann(PDEData3d &myPDE,EulerData3D *myEuler)
{
    
    if (myRiemannMap.find(myRiemannName) != myRiemannMap.end())
        RiemannSolver(myRiemannName,
                      myEuler->getJac() ,
                      myEuler->getXi_xyz(),
                      myEuler->getEta_xyz(),
                      myEuler->getZeta_xyz(),
                      myEuler->getPrimatives(),
                      myEuler->getAcoustics(),
                      myPDE.getQ(),
                      myEuler->getFluxXi(),
                      myEuler->getFluxEta(),
                      myEuler->getFluxZeta());
    else
        std::cout << "Riemann Solver " << myRiemannName << " does not exist!" << std::endl;
}


void nuc3d::physicsModel::getMiu(Field &T,
                                 Field &miu,
                                 Field &coeff)

{
    (this->*myVisModelMap[myVisModelName])(T,
                                           miu,
                                           coeff,
                                           myModelParameters["Reynolds"],
                                           myModelParameters["Mach"],
                                           myModelParameters["Pt"],
                                           myModelParameters["Gamma"],
                                           myModelParameters["T_ref"],
                                           myModelParameters["T_inf"]);
}

void nuc3d::physicsModel::sutherland(const Field &T,
                                     Field &miu,
                                     Field &coeff,
                                     double Reynolds,
                                     double Mach,
                                     double pt,
                                     double gamma,
                                     double T_ref,
                                     double T_inf)
{
    int nx=T.getSizeX();
    int ny=T.getSizeY();
    int nz=T.getSizeZ();
    
    double *pT=T.getDataPtr();
    double *pMiu=miu.getDataPtr();
    double *pCoeff=coeff.getDataPtr();
    
    for (int k = 0; k < nz; k++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                int idx=nx*ny*k+nx*j+i;
                double T_local=pT[idx];
                double non_dim_T_ref=T_ref/T_inf;
                double miu0=(1.0+non_dim_T_ref)/(T_local+non_dim_T_ref)*pow(T_local,1.5)/Reynolds;
                double coeff0=miu0/((gamma-1.0)*Mach*Mach*pt);
                
                pMiu[idx]=miu0;
                pCoeff[idx]=coeff0;
            }
        }
    }
}

void nuc3d::physicsModel::constant(const Field &T,
                                   Field &miu,
                                   Field &coeff,
                                   double Reynolds,
                                   double Mach,
                                   double pt,
                                   double gamma,
                                   double T_ref,
                                   double T_inf)
{
    
}

void nuc3d::physicsModel::initial(PDEData3d &myPDE,std::shared_ptr<EulerData3D> myEuler)
{
    prim2con(myEoSName,
             myEuler->jacobian,
             myEuler->W_Euler,
             myPDE.getQ());
}


void nuc3d::physicsModel::con2prim(const std::string &EoSName,
                                   const Field &Jacobian,
                                   VectorField &Q_vec,
                                   VectorField &W_vec,
                                   VectorField &W0_vec
                                   )
{
    double rho, u, v, w, E, p;
    double T, e, alpha;
    double jac;
    
    int nx = Jacobian.getSizeX();
    int ny = Jacobian.getSizeY();
    int nz = Jacobian.getSizeZ();
    
    auto iterRho = Q_vec.begin();
    auto iterRhoU = Q_vec.begin() + 1;
    auto iterRhoV = Q_vec.begin() + 2;
    auto iterRhoW = Q_vec.begin() + 3;
    auto iterRhoE = Q_vec.begin() + 4;
    
    auto iterRho0 = W_vec.begin();
    auto iterU0   = W_vec.begin() + 1;
    auto iterV0   = W_vec.begin() + 2;
    auto iterW0   = W_vec.begin() + 3;
    auto iterP0   = W_vec.begin() + 4;
    
    auto iterT0 = W0_vec.begin();
    auto iterE0   = W0_vec.begin() + 1;
    auto iterAlpha0   = W0_vec.begin() + 2;
    
    double *pRho=iterRho->getDataPtr();
    double *pRhoU=iterRhoU->getDataPtr();
    double *pRhoV=iterRhoV->getDataPtr();
    double *pRhoW=iterRhoW->getDataPtr();
    double *pRhoE=iterRhoE->getDataPtr();
    double *pjac=Jacobian.getDataPtr();
    
    double *pRho0=iterRho0->getDataPtr();
    double *pU0=iterU0->getDataPtr();
    double *pV0=iterV0->getDataPtr();
    double *pW0=iterW0->getDataPtr();
    double *pP0=iterP0->getDataPtr();
    
    double *pT0=iterT0->getDataPtr();
    double *pE0=iterE0->getDataPtr();
    double *pAlpha0=iterAlpha0->getDataPtr();
    double gamma = myModelParameters["Gamma"];
    for (int k = 0; k < nz; ++k)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                int idx=nx*ny*k+nx*j+i;
                
                jac=pjac[idx];
                rho = pRho[idx]*jac;
                u = pRhoU[idx]/ rho*jac;
                v = pRhoV[idx] / rho*jac;
                w = pRhoW[idx] / rho*jac;
                E = pRhoE[idx]*jac;
                
                double e0=E-0.5*rho*(u*u+v*v+w*w);
                
                if(e0<0.0)
                {
                    //                    if((i!=0)&&(i!=(nx-1)))
//                    {
//                        E= (pRhoE[idx-1]/pjac[idx-1]
//                            +pRhoE[idx+1]/pjac[idx+1])*0.5;
//                        
//                        rho=( pRho[idx-1]/pjac[idx-1]
//                             +pRho[idx+1]/pjac[idx+1])*0.5;
//                        
//                        u= (pRhoU[idx-1]/pjac[idx-1]
//                            +pRhoU[idx+1]/pjac[idx+1])*0.5/ rho;
//                        
//                        v= (pRhoV[idx-1]/pjac[idx-1]
//                            +pRhoV[idx+1]/pjac[idx+1])*0.5/ rho;
//                        
//                        w= (pRhoW[idx-1]/pjac[idx-1]
//                            +pRhoW[idx+1]/pjac[idx+1])*0.5/ rho;
//                    }
//                    else if((j!=0)&&(j!=(ny-1)))
//                    {
//                        E= (pRhoE[idx-nx]/pjac[idx-nx]
//                            +pRhoE[idx+nx]/pjac[idx+nx])*0.5;
//                        
//                        rho=( pRho[idx-nx]/pjac[idx-nx]
//                             +pRho[idx+nx]/pjac[idx+nx])*0.5;
//                        
//                        u= (pRhoU[idx-nx]/pjac[idx-nx]
//                            +pRhoU[idx+nx]/pjac[idx+nx])*0.5/ rho;
//                        
//                        v= (pRhoV[idx-nx]/pjac[idx-nx]
//                            +pRhoV[idx+nx]/pjac[idx+nx])*0.5/ rho;
//                        
//                        w= (pRhoW[idx-nx]/pjac[idx-nx]
//                            +pRhoW[idx+nx]/pjac[idx+nx])*0.5/ rho;}
//                    else if((k!=0)&&(k!=(nz-1)))
//                    {
//                        E= (pRhoE[idx-nx*ny]/pjac[idx-nx*ny]
//                            +pRhoE[idx+nx*ny]/pjac[idx+nx*ny])*0.5;
//                        
//                        rho=( pRho[idx-nx*ny]/pjac[idx-nx*ny]
//                             +pRho[idx+nx*ny]/pjac[idx+nx*ny])*0.5;
//                        
//                        u= (pRhoU[idx-nx*ny]/pjac[idx-nx*ny]
//                            +pRhoU[idx+nx*ny]/pjac[idx+nx*ny])*0.5/ rho;
//                        
//                        v= (pRhoV[idx-nx*ny]/pjac[idx-nx*ny]
//                            +pRhoV[idx+nx*ny]/pjac[idx+nx*ny])*0.5/ rho;
//                        
//                        w= (pRhoW[idx-nx*ny]/pjac[idx-nx*ny]
//                            +pRhoW[idx+nx*ny]/pjac[idx+nx*ny])*0.5/ rho;}
//                    else
//                    {
                        E=(pow(rho,gamma)/(gamma-1.0)+0.5*rho*(u*u+v*v+w*w));
                    //}

                    pRhoE[idx]= E/jac;
                }
                
                (this->*myEosFWDMap[EoSName])(rho,
                                              u,
                                              v,
                                              w,
                                              E,
                                              p,
                                              T,
                                              e,
                                              alpha);
                
                
                pRho0[idx]=rho;
                pU0[idx]=u;
                pV0[idx]=v;
                pW0[idx]=w;
                pP0[idx]=p;
                
                pT0[idx]=T;
                pE0[idx]=e;
                pAlpha0[idx]=alpha;
            }
        }
    }
}

void nuc3d::physicsModel::prim2con(const std::string &EosName,
                                   const Field &Jacobian,
                                   const VectorField &W_vec,
                                   VectorField &Q_vec)
{
    double rho, u, v, w, p;
    double rho0,rhou0,rhov0,rhow0,E;
    double jac;
    
    int nx = Jacobian.getSizeX();
    int ny = Jacobian.getSizeY();
    int nz = Jacobian.getSizeZ();
    
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                auto iterRho = W_vec.begin();
                auto iterU = W_vec.begin() + 1;
                auto iterV = W_vec.begin() + 2;
                auto iterW = W_vec.begin() + 3;
                auto iterP = W_vec.begin() + 4;
                
                rho = iterRho->getValue(i, j, k);
                u  = iterU->getValue(i, j, k);
                v  = iterV->getValue(i, j, k);
                w  = iterW->getValue(i, j, k);
                p  = iterP->getValue(i, j, k);
                jac= Jacobian.getValue(i,j,k);
                
                (this->*myEosBWDMap[myEoSName])(
                                                rho,
                                                u,
                                                v,
                                                w,
                                                p,
                                                rho0,
                                                rhou0,
                                                rhov0,
                                                rhow0,
                                                E);
                
                
                auto iterRho0 = Q_vec.begin();
                auto iterU0 = Q_vec.begin() + 1;
                auto iterV0 = Q_vec.begin() + 2;
                auto iterW0 = Q_vec.begin() + 3;
                auto iterE0 = Q_vec.begin() + 4;
                
                iterRho0->setValue(i, j, k, rho0/jac);
                iterU0->setValue(i, j, k, rhou0/jac);
                iterV0->setValue(i, j, k, rhov0/jac);
                iterW0->setValue(i, j, k, rhow0/jac);
                iterE0->setValue(i, j, k, E/jac);
                
                for (auto iter = iterE0 + 1; iter != Q_vec.end(); ++iter)
                    iter->setValue(i, j, k, Q_vec[iter - Q_vec.end()].getValue(i, j, k) * rho/jac);
            }
    
}

void nuc3d::physicsModel::RiemannSolver(const std::string &SolverName,
                                        const Field &Jacobian,
                                        const VectorField &xi_xyz,
                                        const VectorField &eta_xyz,
                                        const VectorField &zeta_xyz,
                                        const VectorField &W_vec,
                                        const VectorField &W0_vec,
                                        const VectorField &Q_vec,
                                        EulerFlux& myFlux_xi,
                                        EulerFlux& myFlux_eta,
                                        EulerFlux& myFlux_zeta)
{
    int nx = Jacobian.getSizeX();
    int ny = Jacobian.getSizeY();
    int nz = Jacobian.getSizeZ();
    
    double *Rho=Q_vec[0].getDataPtr();
    double *RhoU=Q_vec[1].getDataPtr();
    double *RhoV=Q_vec[2].getDataPtr();
    double *RhoW=Q_vec[3].getDataPtr();
    double *RhoE=Q_vec[4].getDataPtr();
    
    double *rho=W_vec[0].getDataPtr();
    double *u=W_vec[1].getDataPtr();
    double *v=W_vec[2].getDataPtr();
    double *w=W_vec[3].getDataPtr();
    double *p=W_vec[4].getDataPtr();
    
    double *T=W0_vec[0].getDataPtr();
    double *e=W0_vec[1].getDataPtr();
    double *alpha=W0_vec[2].getDataPtr();
    
    double *jac=Jacobian.getDataPtr();
    double *xi_x=xi_xyz[0].getDataPtr();
    double *xi_y=xi_xyz[1].getDataPtr();
    double *xi_z=xi_xyz[2].getDataPtr();
    
    double *eta_x=eta_xyz[0].getDataPtr();
    double *eta_y=eta_xyz[1].getDataPtr();
    double *eta_z=eta_xyz[2].getDataPtr();
    
    double *zeta_x=zeta_xyz[0].getDataPtr();
    double *zeta_y=zeta_xyz[1].getDataPtr();
    double *zeta_z=zeta_xyz[2].getDataPtr();
    
    double U0,V0,W0;
    double theta_xi,theta_eta,theta_zeta;
    double alpha_xi,alpha_eta,alpha_zeta;
    
    double jac0;
    double xi_x0,xi_y0,xi_z0;
    double eta_x0,eta_y0,eta_z0;
    double zeta_x0,zeta_y0,zeta_z0;
    
    double Rho0,Rhou0,Rhov0,Rhow0,Rhoe0;
    double rho0,u0,v0,w0,p0,T0,e0,alpha0;
    
    double MaxEigen_xi=0.0;
    double MaxEigen_eta=0.0;
    double MaxEigen_zeta=0.0;
    double fluxp[5], fluxn[5];
    
    double *flux_xi_l[5];
    double *flux_xi_r[5];
    
    double *flux_eta_l[5];
    double *flux_eta_r[5];
    
    double *flux_zeta_l[5];
    double *flux_zeta_r[5];
    
    flux_xi_l[0]=myFlux_xi.FluxL[0].getDataPtr();
    flux_xi_l[1]=myFlux_xi.FluxL[1].getDataPtr();
    flux_xi_l[2]=myFlux_xi.FluxL[2].getDataPtr();
    flux_xi_l[3]=myFlux_xi.FluxL[3].getDataPtr();
    flux_xi_l[4]=myFlux_xi.FluxL[4].getDataPtr();
    
    flux_xi_r[0]=myFlux_xi.FluxR[0].getDataPtr();
    flux_xi_r[1]=myFlux_xi.FluxR[1].getDataPtr();
    flux_xi_r[2]=myFlux_xi.FluxR[2].getDataPtr();
    flux_xi_r[3]=myFlux_xi.FluxR[3].getDataPtr();
    flux_xi_r[4]=myFlux_xi.FluxR[4].getDataPtr();
    
    flux_eta_l[0]=myFlux_eta.FluxL[0].getDataPtr();
    flux_eta_l[1]=myFlux_eta.FluxL[1].getDataPtr();
    flux_eta_l[2]=myFlux_eta.FluxL[2].getDataPtr();
    flux_eta_l[3]=myFlux_eta.FluxL[3].getDataPtr();
    flux_eta_l[4]=myFlux_eta.FluxL[4].getDataPtr();
    
    flux_eta_r[0]=myFlux_eta.FluxR[0].getDataPtr();
    flux_eta_r[1]=myFlux_eta.FluxR[1].getDataPtr();
    flux_eta_r[2]=myFlux_eta.FluxR[2].getDataPtr();
    flux_eta_r[3]=myFlux_eta.FluxR[3].getDataPtr();
    flux_eta_r[4]=myFlux_eta.FluxR[4].getDataPtr();
    
    flux_zeta_l[0]=myFlux_zeta.FluxL[0].getDataPtr();
    flux_zeta_l[1]=myFlux_zeta.FluxL[1].getDataPtr();
    flux_zeta_l[2]=myFlux_zeta.FluxL[2].getDataPtr();
    flux_zeta_l[3]=myFlux_zeta.FluxL[3].getDataPtr();
    flux_zeta_l[4]=myFlux_zeta.FluxL[4].getDataPtr();
    
    flux_zeta_r[0]=myFlux_zeta.FluxR[0].getDataPtr();
    flux_zeta_r[1]=myFlux_zeta.FluxR[1].getDataPtr();
    flux_zeta_r[2]=myFlux_zeta.FluxR[2].getDataPtr();
    flux_zeta_r[3]=myFlux_zeta.FluxR[3].getDataPtr();
    flux_zeta_r[4]=myFlux_zeta.FluxR[4].getDataPtr();
    
    
    for (int k = 0; k < nz; ++k)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                
                jac0=jac[idx_xi];
                
                xi_x0=xi_x[idx_xi];
                xi_y0=xi_y[idx_xi];
                xi_z0=xi_z[idx_xi];
                
                rho0=rho[idx_xi];
                u0=u[idx_xi];
                v0=v[idx_xi];
                w0=w[idx_xi];
                p0=p[idx_xi];
                T0=T[idx_xi];
                e0=e[idx_xi];
                alpha0=alpha[idx_xi];
                
                Rho0=Rho[idx_xi];
                Rhou0=RhoU[idx_xi];
                Rhov0=RhoV[idx_xi];
                Rhow0=RhoW[idx_xi];
                Rhoe0=RhoE[idx_xi];
                
                theta_xi=sqrt(xi_x0*xi_x0
                              + xi_y0*xi_y0
                              + xi_z0*xi_z0);
                
                U0=(xi_x0*u0 + xi_y0*v0 + xi_z0*w0);
                
                alpha0 = alpha[idx_xi];
                
                alpha_xi = alpha0*theta_xi;
                (this->*myRiemannMap[SolverName])(U0,
                                                  alpha_xi,
                                                  Rho0,
                                                  Rhou0,
                                                  Rhov0,
                                                  Rhow0,
                                                  Rhoe0,
                                                  rho0,
                                                  u0,
                                                  v0,
                                                  w0,
                                                  p0,
                                                  xi_x0,
                                                  xi_y0,
                                                  xi_z0,
                                                  jac0,
                                                  fluxp,
                                                  fluxn);
                
                flux_xi_l[0][idx_xi]=fluxp[0];
                flux_xi_l[1][idx_xi]=fluxp[1];
                flux_xi_l[2][idx_xi]=fluxp[2];
                flux_xi_l[3][idx_xi]=fluxp[3];
                flux_xi_l[4][idx_xi]=fluxp[4];
                
                flux_xi_r[0][idx_xi]=fluxn[0];
                flux_xi_r[1][idx_xi]=fluxn[1];
                flux_xi_r[2][idx_xi]=fluxn[2];
                flux_xi_r[3][idx_xi]=fluxn[3];
                flux_xi_r[4][idx_xi]=fluxn[4];
                
                MaxEigen_xi=std::max(std::abs(U0)+alpha_xi ,MaxEigen_xi);
                
            }
        }
    }
    
    for (int k = 0; k < nz; ++k)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                int idx_eta=ny*nz*i+ny*k+j;
                
                jac0=jac[idx_xi];
                eta_x0=eta_x[idx_xi];
                eta_y0=eta_y[idx_xi];
                eta_z0=eta_z[idx_xi];
                
                rho0=rho[idx_xi];
                u0=u[idx_xi];
                v0=v[idx_xi];
                w0=w[idx_xi];
                p0=p[idx_xi];
                T0=T[idx_xi];
                e0=e[idx_xi];
                alpha0=alpha[idx_xi];
                
                Rho0=Rho[idx_xi];
                Rhou0=RhoU[idx_xi];
                Rhov0=RhoV[idx_xi];
                Rhow0=RhoW[idx_xi];
                Rhoe0=RhoE[idx_xi];
                
                theta_eta=std::sqrt(eta_x0*eta_x0
                                    + eta_y0*eta_y0
                                    + eta_z0*eta_z0);
                
                V0=(eta_x0*u0 + eta_y0*v0 + eta_z0*w0);
                
                alpha0 = alpha[idx_xi];
                
                alpha_eta = alpha0*theta_eta;
                
                (this->*myRiemannMap[SolverName])(V0,
                                                  alpha_eta,
                                                  Rho0,
                                                  Rhou0,
                                                  Rhov0,
                                                  Rhow0,
                                                  Rhoe0,
                                                  rho0,
                                                  u0,
                                                  v0,
                                                  w0,
                                                  p0,
                                                  eta_x0,
                                                  eta_y0,
                                                  eta_z0,
                                                  jac0,
                                                  fluxp,
                                                  fluxn);
                flux_eta_l[0][idx_eta]=fluxp[0];
                flux_eta_l[1][idx_eta]=fluxp[1];
                flux_eta_l[2][idx_eta]=fluxp[2];
                flux_eta_l[3][idx_eta]=fluxp[3];
                flux_eta_l[4][idx_eta]=fluxp[4];
                
                flux_eta_r[0][idx_eta]=fluxn[0];
                flux_eta_r[1][idx_eta]=fluxn[1];
                flux_eta_r[2][idx_eta]=fluxn[2];
                flux_eta_r[3][idx_eta]=fluxn[3];
                flux_eta_r[4][idx_eta]=fluxn[4];
                
                
                MaxEigen_eta=std::max(std::abs(V0)+alpha_eta,MaxEigen_eta);
                
            }
        }
    }
    
    for (int k = 0; k < nz; ++k)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                int idx_xi=nx*ny*k+nx*j+i;
                int idx_zeta=nz*nx*j+nz*i+k;
                
                jac0=jac[idx_xi];
                
                zeta_x0=zeta_x[idx_xi];
                zeta_y0=zeta_y[idx_xi];
                zeta_z0=zeta_z[idx_xi];
                
                rho0=rho[idx_xi];
                u0=u[idx_xi];
                v0=v[idx_xi];
                w0=w[idx_xi];
                p0=p[idx_xi];
                T0=T[idx_xi];
                e0=e[idx_xi];
                alpha0=alpha[idx_xi];
                
                Rho0=Rho[idx_xi];
                Rhou0=RhoU[idx_xi];
                Rhov0=RhoV[idx_xi];
                Rhow0=RhoW[idx_xi];
                Rhoe0=RhoE[idx_xi];
                
                theta_zeta=sqrt(zeta_x0*zeta_x0
                                + zeta_y0*zeta_y0
                                + zeta_z0*zeta_z0);
                
                
                W0=(zeta_x0*u0 + zeta_y0*v0 + zeta_z0*w0);
                
                alpha0 = alpha[idx_xi];
                alpha_zeta = alpha0*theta_zeta;
                
                (this->*myRiemannMap[SolverName])(W0,
                                                  alpha_zeta,
                                                  Rho0,
                                                  Rhou0,
                                                  Rhov0,
                                                  Rhow0,
                                                  Rhoe0,
                                                  rho0,
                                                  u0,
                                                  v0,
                                                  w0,
                                                  p0,
                                                  zeta_x0,
                                                  zeta_y0,
                                                  zeta_z0,
                                                  jac0,
                                                  fluxp,
                                                  fluxn);
                
                flux_zeta_l[0][idx_zeta]=fluxp[0];
                flux_zeta_l[1][idx_zeta]=fluxp[1];
                flux_zeta_l[2][idx_zeta]=fluxp[2];
                flux_zeta_l[3][idx_zeta]=fluxp[3];
                flux_zeta_l[4][idx_zeta]=fluxp[4];
                
                flux_zeta_r[0][idx_zeta]=fluxn[0];
                flux_zeta_r[1][idx_zeta]=fluxn[1];
                flux_zeta_r[2][idx_zeta]=fluxn[2];
                flux_zeta_r[3][idx_zeta]=fluxn[3];
                flux_zeta_r[4][idx_zeta]=fluxn[4];
                
                MaxEigen_zeta=std::max(std::abs(W0)+alpha_zeta,MaxEigen_zeta);
                
            }
        }
    }
    
    myFlux_xi.maxEigen=MaxEigen_xi;
    myFlux_eta.maxEigen=MaxEigen_eta;
    myFlux_zeta.maxEigen=MaxEigen_zeta;
    
};

void nuc3d::physicsModel::RiemannAUSM(const double &U0,
                                      const double &alpha0,
                                      const double &Rho,
                                      const double &RhoU,
                                      const double &RhoV,
                                      const double &RhoW,
                                      const double &RhoE,
                                      const double &rho,
                                      const double &u,
                                      const double &v,
                                      const double &w,
                                      const double &p,
                                      const double &xx_x,
                                      const double &xx_y,
                                      const double &xx_z,
                                      const double &jac,
                                      double *fluxp,
                                      double *fluxn)
{
    double mach=U0/alpha0;
    double machp = getMachL(mach);
    double machn = getMachR(mach);
    
    double p_p = getPressureL(mach, p);
    double p_n = getPressureR(mach, p);
    
    fluxp[0] = machp*alpha0*Rho;
    fluxp[1] = machp*alpha0*RhoU + xx_x*p_p/jac;
    fluxp[2] = machp*alpha0*RhoV + xx_y*p_p/jac;
    fluxp[3] = machp*alpha0*RhoW + xx_z*p_p/jac;
    fluxp[4] = machp*alpha0*RhoE + machp*alpha0*p/jac;
    
    fluxn[0] = machn*alpha0*Rho;
    fluxn[1] = machn*alpha0*RhoU + xx_x*p_n/jac;
    fluxn[2] = machn*alpha0*RhoV + xx_y*p_n/jac;
    fluxn[3] = machn*alpha0*RhoW + xx_z*p_n/jac;
    fluxn[4] = machn*alpha0*RhoE + machn*alpha0*p/jac;
}

double nuc3d::physicsModel::getMachL(const double &mach)
{
    double MachL;
    
    if (std::abs(mach) < 1.0)
        MachL = 0.25*pow((mach + 1.0), 2);
    else
        MachL = 0.50*(mach + std::abs(mach));
    
    return MachL;
    
}

double nuc3d::physicsModel::getMachR(const double &mach)
{
    
    double MachR;
    
    if (std::abs(mach) < 1.0)
        MachR = -0.25*pow((mach - 1.0), 2);
    else
        MachR = 0.50*(mach - std::abs(mach));
    
    return MachR;
    
}

double nuc3d::physicsModel::getPressureL(const double &mach, const double &p)
{
    double pressureL;
    if (std::abs(mach) < 1.0)
        pressureL = p*(0.25*pow((mach + 1.0), 2)*(2.0 - mach));
    else
        pressureL = 0.50*p*(mach + std::abs(mach)) / mach;
    
    return pressureL;
}

double nuc3d::physicsModel::getPressureR(const double &mach, const double &p)
{
    double pressureR;
    if (std::abs(mach) < 1.0)
        pressureR = p*(0.25*pow(mach - 1.0, 2)*(2.0 + mach));
    else
        pressureR = 0.5*p*(mach - std::abs(mach)) / mach;
    
    return pressureR;
}

void nuc3d::physicsModel::RiemannAUSMp(const double &U0,
                                       const double &alpha0,
                                       const double &Rho,
                                       const double &RhoU,
                                       const double &RhoV,
                                       const double &RhoW,
                                       const double &RhoE,
                                       const double &rho,
                                       const double &u,
                                       const double &v,
                                       const double &w,
                                       const double &p,
                                       const double &xx_x,
                                       const double &xx_y,
                                       const double &xx_z,
                                       const double &jac,
                                       double *fluxp,
                                       double *fluxn)
{
    double mach=U0/alpha0;
    double machp = getMachLp(mach);
    double machn = getMachRp(mach);
    
    double p_p = getPressureLp(mach, p);
    double p_n = getPressureRp(mach, p);
    
    fluxp[0] = machp*alpha0*Rho;
    fluxp[1] = machp*alpha0*RhoU + xx_x*p_p/jac;
    fluxp[2] = machp*alpha0*RhoV + xx_y*p_p/jac;
    fluxp[3] = machp*alpha0*RhoW + xx_z*p_p/jac;
    fluxp[4] = machp*alpha0*RhoE + machp*alpha0*p/jac;
    
    fluxn[0] = machn*alpha0*Rho;
    fluxn[1] = machn*alpha0*RhoU + xx_x*p_n/jac;
    fluxn[2] = machn*alpha0*RhoV + xx_y*p_n/jac;
    fluxn[3] = machn*alpha0*RhoW + xx_z*p_n/jac;
    fluxn[4] = machn*alpha0*RhoE + machn*alpha0*p/jac;
}

double nuc3d::physicsModel::getMachLp(const double &mach)
{
    double MachL;
    
    if (std::abs(mach) < 1.0)
        MachL = 0.25*pow((mach + 1.0), 2)+0.125*pow((mach*mach-1), 2);
    else
        MachL = 0.50*(mach + std::abs(mach));
    
    return MachL;
    
}

double nuc3d::physicsModel::getMachRp(const double &mach)
{
    
    double MachR;
    
    if (std::abs(mach) < 1.0)
        MachR = -0.25*pow((mach - 1.0), 2)-0.125*pow((mach*mach-1), 2);
    else
        MachR = 0.50*(mach - std::abs(mach));
    
    return MachR;
    
}

double nuc3d::physicsModel::getPressureLp(const double &mach, const double &p)
{
    double pressureL;
    if (std::abs(mach) < 1.0)
        pressureL = p*(0.25*pow((mach + 1.0), 2)*(2.0 - mach)+0.1875*mach*pow((mach*mach-1.0),2));
    else
        pressureL = 0.50*p*(mach + std::abs(mach)) / mach;
    
    return pressureL;
}

double nuc3d::physicsModel::getPressureRp(const double &mach, const double &p)
{
    double pressureR;
    if (std::abs(mach) < 1.0)
        pressureR = p*(0.25*pow(mach - 1.0, 2)*(2.0 + mach)-0.1875*mach*pow((mach*mach-1.0),2));
    else
        pressureR = 0.5*p*(mach - std::abs(mach)) / mach;
    
    return pressureR;
}

void nuc3d::physicsModel::RiemannLF(const double &U0,
                                    const double &alpha0,
                                    const double &Rho,
                                    const double &RhoU,
                                    const double &RhoV,
                                    const double &RhoW,
                                    const double &RhoE,
                                    const double &rho,
                                    const double &u,
                                    const double &v,
                                    const double &w,
                                    const double &p,
                                    const double &xx_x,
                                    const double &xx_y,
                                    const double &xx_z,
                                    const double &jac,
                                    double *fluxp,
                                    double *fluxn)
{
    double flux[5];
    double MaxEigenLocal=std::abs(U0)+alpha0;
    
    flux[0]=U0*Rho;
    flux[1]=U0*RhoU+ xx_x*p/jac;
    flux[2]=U0*RhoV+ xx_y*p/jac;
    flux[3]=U0*RhoW+ xx_z*p/jac;
    flux[4]=U0*(RhoE+p/jac);
    
    fluxp[0] = 0.5*(flux[0]+MaxEigenLocal*Rho);
    fluxp[1] = 0.5*(flux[1]+MaxEigenLocal*RhoU);
    fluxp[2] = 0.5*(flux[2]+MaxEigenLocal*RhoV);
    fluxp[3] = 0.5*(flux[3]+MaxEigenLocal*RhoW);
    fluxp[4] = 0.5*(flux[4]+MaxEigenLocal*RhoE);
    
    fluxn[0] = 0.5*(flux[0]-MaxEigenLocal*Rho);
    fluxn[1] = 0.5*(flux[1]-MaxEigenLocal*RhoU);
    fluxn[2] = 0.5*(flux[2]-MaxEigenLocal*RhoV);
    fluxn[3] = 0.5*(flux[3]-MaxEigenLocal*RhoW);
    fluxn[4] = 0.5*(flux[4]-MaxEigenLocal*RhoE);
}

void nuc3d::physicsModel::EoSIdealGasFWD(const double &rho,
                                         const double &u,
                                         const double &v,
                                         const double &w,
                                         const double &E,
                                         double &p,
                                         double &T,
                                         double &e,
                                         double &alpha)
{
    double gamma = myModelParameters["Gamma"];
    double Mach = myModelParameters["Mach"];
    
    p = EoSIdealGasgetPressure(rho, u, v, w, E, Mach, gamma);
    T = EoSIdealGasgetTemp(rho, u, v, w, E,p, Mach, gamma);
    e = EoSIdealGasgetE0(rho, u, v, w, E,p, Mach, gamma);
    alpha =EoSIdealGasgetAlpha(rho, u, v, w, E, p,Mach, gamma);
};

double nuc3d::physicsModel::EoSIdealGasgetPressure(const double &rho,
                                                   const double &u,
                                                   const double &v,
                                                   const double &w,
                                                   const double &E,
                                                   const double &mach,
                                                   const double &gamma)
{
    double temp=(E - 0.5*rho*(u*u + v*v + w*w))*(gamma - 1.0);
    if(temp<0.0)
    {
        int myID;
        MPI_Comm_rank(MPI_COMM_WORLD, &myID);
        std::cout<<"negative pressure!!!"<<" Proc ID = "<<myID<<std::endl;
        exit(-1);
        //temp=std::pow(rho,gamma);
    }
    return (E - 0.5*rho*(u*u + v*v + w*w))*(gamma - 1.0);
}

double nuc3d::physicsModel::EoSIdealGasgetTemp(const double &rho,
                                               const double &u,
                                               const double &v,
                                               const double &w,
                                               const double &E,
                                               const double &p,
                                               const double &mach,
                                               const double &gamma)
{
    return gamma*mach*mach*p / rho;
}

double nuc3d::physicsModel::EoSIdealGasgetE0(const double &rho,
                                             const double &u,
                                             const double &v,
                                             const double &w,
                                             const double &E,
                                             const double &p,
                                             const double &mach,
                                             const double &gamma)
{
    return p / rho / (gamma - 1.0);
}

double nuc3d::physicsModel::EoSIdealGasgetAlpha(const double &rho,
                                                const double &u,
                                                const double &v,
                                                const double &w,
                                                const double &E,
                                                const double &p,
                                                const double &mach,
                                                const double &gamma)
{
    return std::sqrt(std::abs(gamma*p/rho));
}

void nuc3d::physicsModel::EoSIdealGasBWD(
                                         const double &rho,
                                         const double &u,
                                         const double &v,
                                         const double &w,
                                         const double &p,
                                         double & _rho,
                                         double & _rhou,
                                         double & _rhov,
                                         double & _rhow,
                                         double & _E)
{
    double gamma = myModelParameters["Gamma"];
    
    _rho = rho;
    _rhou= rho*u;
    _rhov= rho*v;
    _rhow= rho*w;
    _E   = p/(gamma-1.0)+0.5*rho*(u*u + v*v + w*w);
    
    
};


void nuc3d::physicsModel::EoSJWLFWD(const double &rho,
                                    const double &u,
                                    const double &v,
                                    const double &w,
                                    const double &E, //E=p/(gamma-1)+1/2*rho*(u^2+v^2+w^2)
                                    double &p,
                                    double &T,
                                    double &e,
                                    double &alpha)
{
};
void nuc3d::physicsModel::EoSJWLBWD(
                                    const double &rho,
                                    const double &u,
                                    const double &v,
                                    const double &w,
                                    const double &p,
                                    double & _rho,
                                    double & _rhou,
                                    double & _rhov,
                                    double & _rhow,
                                    double & _E)
{
};

void nuc3d::physicsModel::prim2conPoint(const double &rho,
                                        const double &u,
                                        const double &v,
                                        const double &w,
                                        const double &p,
                                        double &T,
                                        double &e,
                                        double &alpha)
{
    double _rho;
    double _rhou;
    double _rhov;
    double _rhow;
    double _E;
    double _P;
    
    (this->*myEosBWDMap[myEoSName])(rho,
                                    u,
                                    v,
                                    w,
                                    p,
                                    _rho,
                                    _rhou,
                                    _rhov,
                                    _rhow,
                                    _E);
    
    (this->*myEosFWDMap[myEoSName])(rho,
                                    u,
                                    v,
                                    w,
                                    _E,
                                    _P,
                                    T,
                                    e,
                                    alpha);
    
    
}

void nuc3d::physicsModel::solveRiemannPoint(const std::vector<double> &prim,
                                            const double jac,
                                            const double xx_x,
                                            const double xx_y,
                                            const double xx_z,
                                            std::vector<double> &fluxl,
                                            std::vector<double> &fluxr)
{
    if (myRiemannMap.find(myRiemannName) != myRiemannMap.end())
        (this->*myRiemannPointMap[myRiemannName])(prim,jac,xx_x,xx_y,xx_z,fluxl,fluxr);
    else
        std::cout << "Riemann Solver " << myRiemannName << " does not exist!" << std::endl;
}

void nuc3d::physicsModel::solveRiemannPointAUSM(const std::vector<double> &prim,
                                                const double jac,
                                                const double xx_x,
                                                const double xx_y,
                                                const double xx_z,
                                                std::vector<double> &fluxl,
                                                std::vector<double> &fluxr)
{
    const double &rho=prim[0];
    const double &u=prim[1];
    const double &v=prim[2];
    const double &w=prim[3];
    const double &p=prim[4];
    
    double _rho;
    double _rhou;
    double _rhov;
    double _rhow;
    double _E;
    double _p;
    double _T;
    double _e;
    double _alpha;
    
    
    (this->*myEosBWDMap[myEoSName])(rho,
                                    u,
                                    v,
                                    w,
                                    p,
                                    _rho,
                                    _rhou,
                                    _rhov,
                                    _rhow,
                                    _E);
    
    (this->*myEosFWDMap[myEoSName])(rho,
                                    u,
                                    v,
                                    w,
                                    _E,
                                    _p,
                                    _T,
                                    _e,
                                    _alpha);
    
    double mach, machp, machn;
    double p_p, p_n;
    double U0;
    double theta;
    theta=sqrt(xx_x*xx_x + xx_y*xx_y + xx_z*xx_z);
    double alpha = _alpha*theta;
    
    U0=(xx_x*u + xx_y*v + xx_z*w);
    
    mach = U0 / alpha;
    
    
    machp = getMachL(mach);
    machn = getMachR(mach);
    
    p_p = getPressureL(mach, p);
    p_n = getPressureR(mach, p);
    
    
    
    fluxl[0] = machp*alpha*_rho/jac;
    fluxl[1] = machp*alpha*_rhou/jac + xx_x*p_p/jac;
    fluxl[2] = machp*alpha*_rhov/jac + xx_y*p_p/jac;
    fluxl[3] = machp*alpha*_rhow/jac + xx_z*p_p/jac;
    fluxl[4] = machp*alpha*_E/jac + machp*alpha*p/jac;
    
    fluxr[0] = machn*alpha*_rho/jac;
    fluxr[1] = machn*alpha*_rhou/jac + xx_x*p_n/jac;
    fluxr[2] = machn*alpha*_rhov/jac + xx_y*p_n/jac;
    fluxr[3] = machn*alpha*_rhow/jac + xx_z*p_n/jac;
    fluxr[4] = machn*alpha*_E/jac + machn*alpha*p/jac;
    
    
    
}

void nuc3d::physicsModel::solveRiemannPointAUSMp(const std::vector<double> &prim,
                                                 const double jac,
                                                 const double xx_x,
                                                 const double xx_y,
                                                 const double xx_z,
                                                 std::vector<double> &fluxl,
                                                 std::vector<double> &fluxr)
{
    const double &rho=prim[0];
    const double &u=prim[1];
    const double &v=prim[2];
    const double &w=prim[3];
    const double &p=prim[4];
    
    double _rho;
    double _rhou;
    double _rhov;
    double _rhow;
    double _E;
    double _p;
    double _T;
    double _e;
    double _alpha;
    
    
    (this->*myEosBWDMap[myEoSName])(rho,
                                    u,
                                    v,
                                    w,
                                    p,
                                    _rho,
                                    _rhou,
                                    _rhov,
                                    _rhow,
                                    _E);
    
    (this->*myEosFWDMap[myEoSName])(rho,
                                    u,
                                    v,
                                    w,
                                    _E,
                                    _p,
                                    _T,
                                    _e,
                                    _alpha);
    
    double mach, machp, machn;
    double p_p, p_n;
    double U0;
    double theta;
    theta=sqrt(xx_x*xx_x + xx_y*xx_y + xx_z*xx_z);
    double alpha = _alpha*theta;
    
    U0=(xx_x*u + xx_y*v + xx_z*w);
    
    mach = U0 / alpha;
    
    
    machp = getMachLp(mach);
    machn = getMachRp(mach);
    
    p_p = getPressureLp(mach, p);
    p_n = getPressureRp(mach, p);
    
    
    
    fluxl[0] = machp*alpha*_rho/jac;
    fluxl[1] = machp*alpha*_rhou/jac + xx_x*p_p/jac;
    fluxl[2] = machp*alpha*_rhov/jac + xx_y*p_p/jac;
    fluxl[3] = machp*alpha*_rhow/jac + xx_z*p_p/jac;
    fluxl[4] = machp*alpha*_E/jac + machp*alpha*p/jac;
    
    fluxr[0] = machn*alpha*_rho/jac;
    fluxr[1] = machn*alpha*_rhou/jac + xx_x*p_n/jac;
    fluxr[2] = machn*alpha*_rhov/jac + xx_y*p_n/jac;
    fluxr[3] = machn*alpha*_rhow/jac + xx_z*p_n/jac;
    fluxr[4] = machn*alpha*_E/jac + machn*alpha*p/jac;
    
    
    
}

void nuc3d::physicsModel::solveRiemannPointLF(const std::vector<double> &prim,
                                              const double jac,
                                              const double xx_x,
                                              const double xx_y,
                                              const double xx_z,
                                              std::vector<double> &fluxl,
                                              std::vector<double> &fluxr)
{
    const double &rho=prim[0];
    const double &u=prim[1];
    const double &v=prim[2];
    const double &w=prim[3];
    const double &p=prim[4];
    double flux[5];
    
    double _rho;
    double _rhou;
    double _rhov;
    double _rhow;
    double _E;
    double _p;
    double _T;
    double _e;
    double _alpha;
    
    
    (this->*myEosBWDMap[myEoSName])(rho,
                                    u,
                                    v,
                                    w,
                                    p,
                                    _rho,
                                    _rhou,
                                    _rhov,
                                    _rhow,
                                    _E);
    
    (this->*myEosFWDMap[myEoSName])(rho,
                                    u,
                                    v,
                                    w,
                                    _E,
                                    _p,
                                    _T,
                                    _e,
                                    _alpha);
    
    double mach;
    double U0;
    double theta;
    theta=sqrt(xx_x*xx_x + xx_y*xx_y + xx_z*xx_z);
    double alpha = _alpha*theta;
    
    U0=(xx_x*u + xx_y*v + xx_z*w);
    
    mach = U0 / alpha;
    double MaxEigenLocal=std::abs(U0)+alpha;
    
    flux[0]=U0*_rho/jac;
    flux[1]=U0*_rhou/jac+ xx_x*p/jac;
    flux[2]=U0*_rhov/jac+ xx_y*p/jac;
    flux[3]=U0*_rhow/jac+ xx_z*p/jac;
    flux[4]=U0*(_E/jac+p/jac);
    
    fluxl[0] = 0.5*(flux[0]+MaxEigenLocal*_rho/jac);
    fluxl[1] = 0.5*(flux[1]+MaxEigenLocal*_rhou/jac);
    fluxl[2] = 0.5*(flux[2]+MaxEigenLocal*_rhov/jac);
    fluxl[3] = 0.5*(flux[3]+MaxEigenLocal*_rhow/jac);
    fluxl[4] = 0.5*(flux[4]+MaxEigenLocal*_E/jac);
    
    fluxr[0] = 0.5*(flux[0]-MaxEigenLocal*_rho/jac);
    fluxr[1] = 0.5*(flux[1]-MaxEigenLocal*_rhou/jac);
    fluxr[2] = 0.5*(flux[2]-MaxEigenLocal*_rhov/jac);
    fluxr[3] = 0.5*(flux[3]-MaxEigenLocal*_rhow/jac);
    fluxr[4] = 0.5*(flux[4]-MaxEigenLocal*_E/jac);
}
/**************************************************************************************
 End of definition
 
 **************************************************************************************/
#endif
