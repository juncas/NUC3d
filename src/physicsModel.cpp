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
                 { "AUSM",&nuc3d::physicsModel::RiemannAUSM },
                 { "LF",&nuc3d::physicsModel::RiemannLF }
             }
             ),
myModelParameters(
                  {
                      {"Reynolds",1000},
                      {"Mach",0.5},
                      {"Pt",0.72},
                      {"Gamma",1.4},
                      {"T_ref",110.4},
                      {"T_inf",288.0}
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
    RiemannSolver(myRiemannName,
                  myEuler->jacobian,
                  myEuler->xi_xyz,
                  myEuler->W_Euler,
                  myEuler->W0_Euler,
                  myPDE.getQ(),
                  myEuler->Flux_xi);
    
    RiemannSolver(myRiemannName,
                  myEuler->jacobian,
                  myEuler->eta_xyz,
                  myEuler->W_Euler,
                  myEuler->W0_Euler,
                  myPDE.getQ(),
                  myEuler->Flux_eta);
    
    RiemannSolver(myRiemannName,
                  myEuler->jacobian,
                  myEuler->zeta_xyz,
                  myEuler->W_Euler,
                  myEuler->W0_Euler,
                  myPDE.getQ(),
                  myEuler->Flux_zeta);
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
                                           myModelParameters["pt"],
                                           myModelParameters["gamma"],
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
                                   const VectorField &Q_vec,
                                   VectorField &W_vec,
                                   VectorField &W0_vec
                                   )
{
    double rho, u, v, w, E, p;
    double T, e, alpha;
    
    int nx = Jacobian.getSizeX();
    int ny = Jacobian.getSizeY();
    int nz = Jacobian.getSizeZ();
    
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                auto iterRho = Q_vec.begin();
                auto iterRhoU = Q_vec.begin() + 1;
                auto iterRhoV = Q_vec.begin() + 2;
                auto iterRhoW = Q_vec.begin() + 3;
                auto iterRhoE = Q_vec.begin() + 4;
                
                double jac=Jacobian.getValue(i, j, k);
                rho = iterRho->getValue(i, j, k)*jac;
                u = iterRhoU->getValue(i, j, k) / rho*jac;
                v = iterRhoV->getValue(i, j, k) / rho*jac;
                w = iterRhoW->getValue(i, j, k) / rho*jac;
                E = iterRhoE->getValue(i, j, k)*jac;
                
                (this->*myEosFWDMap[EoSName])(
                                              rho,
                                              u,
                                              v,
                                              w,
                                              E,
                                              p,
                                              T,
                                              e,
                                              alpha);
                
                
                auto iterRho0 = W_vec.begin();
                auto iterU0   = W_vec.begin() + 1;
                auto iterV0   = W_vec.begin() + 2;
                auto iterW0   = W_vec.begin() + 3;
                auto iterE0   = W_vec.begin() + 4;
                
                iterRho0->setValue(i, j, k, rho);
                iterU0->setValue(i, j, k, u);
                iterV0->setValue(i, j, k, v);
                iterW0->setValue(i, j, k, w);
                iterE0->setValue(i, j, k, p);
                
                for (auto iter = iterE0 + 1; iter != W_vec.end(); ++iter)
                    iter->setValue(i, j, k, Q_vec[iter - W_vec.end()].getValue(i, j, k) / rho*jac);
                
                W0_vec[0].setValue(i, j, k, T);
                W0_vec[1].setValue(i, j, k, e);
                W0_vec[2].setValue(i, j, k, alpha);
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
                                        const VectorField &xx_xyz,
                                        const VectorField &W_vec,
                                        const VectorField &W0_vec,
                                        const VectorField &Q_vec,
                                        EulerFlux& myFlux)
{
    if (myRiemannMap.find(SolverName) != myRiemannMap.end())
        (this->*myRiemannMap[SolverName])(Jacobian, xx_xyz, W_vec, W0_vec, Q_vec, myFlux.FluxL, myFlux.FluxR, myFlux.maxEigen);
    else
        std::cout << "Riemann Solver " << SolverName << " does not exist!" << std::endl;
};

void nuc3d::physicsModel::RiemannAUSM(
                                      const Field &Jacobian,
                                      const VectorField &xx_xyz,
                                      const VectorField &W_vec,
                                      const VectorField &W0_vec,
                                      const VectorField &Q_vec,
                                      VectorField &FluxL,
                                      VectorField &FluxR,
                                      double &MaxEigen)
{
    double rho, u, v, w;
    double alpha;
    double xx_x, xx_y, xx_z, jac;
    double mach, machp, machn;
    double p, p_p, p_n;
    double U0;
    double theta;
    double MaxEigenLocal;
    double fluxp[5], fluxn[5];
    
    int nx = Jacobian.getSizeX();
    int ny = Jacobian.getSizeY();
    int nz = Jacobian.getSizeZ();
    
    MaxEigen=0.0;
    
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                auto iterX = xx_xyz.begin();
                auto iterY = xx_xyz.begin() + 1;
                auto iterZ = xx_xyz.begin() + 2;
                
                xx_x = iterX->getValue(i, j, k);
                xx_y = iterY->getValue(i, j, k);
                xx_z = iterZ->getValue(i, j, k);
                jac = Jacobian.getValue(i, j, k);
                
                rho = W_vec[0].getValue(i, j, k);
                u = W_vec[1].getValue(i, j, k);
                v = W_vec[2].getValue(i, j, k);
                w = W_vec[3].getValue(i, j, k);
                p = W_vec[4].getValue(i, j, k);
                
                theta=sqrt(xx_x*xx_x + xx_y*xx_y + xx_z*xx_z);
                alpha = W0_vec[2].getValue(i, j, k)*theta;
                
                U0=(xx_x*u + xx_y*v + xx_z*w);
                MaxEigenLocal=std::abs(U0)+alpha;
                mach = U0 / alpha;
                
                MaxEigen=(MaxEigen>=MaxEigenLocal)?MaxEigen:MaxEigenLocal;
                
                machp = getMachL(mach);
                machn = getMachR(mach);
                
                p_p = getPressureL(mach, p);
                p_n = getPressureR(mach, p);
                
                
                auto iterRho = Q_vec.begin();
                auto iterU = Q_vec.begin() + 1;
                auto iterV = Q_vec.begin() + 2;
                auto iterW = Q_vec.begin() + 3;
                auto iterE = Q_vec.begin() + 4;
                
                fluxp[0] = machp*alpha*iterRho->getValue(i, j, k);
                fluxp[1] = machp*alpha*iterU->getValue(i, j, k) + xx_x*p_p/jac;
                fluxp[2] = machp*alpha*iterV->getValue(i, j, k) + xx_y*p_p/jac;
                fluxp[3] = machp*alpha*iterW->getValue(i, j, k) + xx_z*p_p/jac;
                fluxp[4] = machp*alpha*iterE->getValue(i, j, k) + machp*alpha*p/jac;
                
                fluxn[0] = machn*alpha*iterRho->getValue(i, j, k);
                fluxn[1] = machn*alpha*iterU->getValue(i, j, k) + xx_x*p_n/jac;
                fluxn[2] = machn*alpha*iterV->getValue(i, j, k) + xx_y*p_n/jac;
                fluxn[3] = machn*alpha*iterW->getValue(i, j, k) + xx_z*p_n/jac;
                fluxn[4] = machn*alpha*iterE->getValue(i, j, k) + machn*alpha*p/jac;
                
                for (auto iter = FluxL.begin(); (iter - FluxL.begin()) != 5; iter++)
                    iter->setValue(i, j, k, fluxp[iter - FluxL.begin()]);
                for (auto iter = FluxL.begin() + 5; iter != FluxL.end(); iter++)
                    iter->setValue(i, j, k, machp*alpha*Q_vec[iter - FluxL.begin()].getValue(i, j, k));
                
                for (auto iter = FluxR.begin(); (iter - FluxR.begin()) != 5; iter++)
                    iter->setValue(i, j, k, fluxn[iter - FluxR.begin()]);
                for (auto iter = FluxR.begin() + 5; iter != FluxR.end(); iter++)
                    iter->setValue(i, j, k, machn*alpha*Q_vec[iter - FluxR.begin()].getValue(i, j, k));
                
            }
};

double nuc3d::physicsModel::getMachL(const double &mach)
{
    double MachL;
    
    if (std::abs(mach) < 1.0)
        MachL = 0.25*pow((mach + 1.0), 2) + 0.125*pow(pow(mach, 2) - 1.0, 2);
    else
        MachL = 0.50*(mach + std::abs(mach));
    
    return MachL;
    
}

double nuc3d::physicsModel::getMachR(const double &mach)
{
    
    double MachR;
    
    if (std::abs(mach) < 1.0)
        MachR = -0.25*pow((mach - 1.0), 2) - 0.125*pow(pow(mach, 2) - 1.0, 2);
    else
        MachR = 0.50*(mach - std::abs(mach));
    
    return MachR;
    
}

double nuc3d::physicsModel::getPressureL(const double &mach, const double &p)
{    double pressureL;
    if (std::abs(mach) < 1.0)
        pressureL = p*(0.25*pow((mach + 1.0), 2)*(2.0 - mach) + 0.1875*mach*pow(pow(mach, 2.0) - 1.0, 2));
    else
        pressureL = 0.50*p*(mach + std::abs(mach)) / mach;
    
    return pressureL;
}

double nuc3d::physicsModel::getPressureR(const double &mach, const double &p)
{
    double pressureR;
    if (std::abs(mach) < 1.0)
        pressureR = p*(0.25*pow(mach - 1.0, 2)*(2.0 + mach) - 0.1875*mach*pow(pow(mach, 2.0) - 1.0, 2));
    else
        pressureR = 0.5*p*(mach - std::abs(mach)) / mach;
    
    return pressureR;
}

void nuc3d::physicsModel::RiemannLF(
                                    const Field &Jacobian,
                                    const VectorField &xx_xyz,
                                    const VectorField &W_vec,
                                    const VectorField &W0_vec,
                                    const VectorField &Q_vec,
                                    VectorField &FluxL,
                                    VectorField &FluxR,
                                    double &MaxEigen)
{
    
};

void nuc3d::physicsModel::EoSIdealGasFWD(const double &rho,
                                         const double &u,
                                         const double &v,
                                         const double &w,
                                          double &E,
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
/**************************************************************************************
 End of definition
 
 **************************************************************************************/
#endif
