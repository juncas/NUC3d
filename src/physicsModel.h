#ifndef physicsModel_h
#define physicsModel_h

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <cmath>
#include "field.h"
#include "euler3d.h"

namespace nuc3d
{
    class PDEData3d;
    class EulerFlux;
    class EulerData3D;
    class EulerReactiveData3D;
    class NaiverStokesData3d;
    class NaiverStokesReactiveData3d;
    
    class physicsModel
    {
        //rho rho*u rho*v rho*w E -> rho u v w p
        typedef void (physicsModel::*pEoSFWD)(const double &,
                                              const double &,
                                              const double &,
                                              const double &,
                                              const double &,
                                              double &,
                                              double &,
                                              double &,
                                              double &);
        
        //rho u v w p -> rho rho*u rho*v rho*w E
        typedef void (physicsModel::*pEoSBWD)(
                                              const double &,
                                              const double &,
                                              const double &,
                                              const double &,
                                              const double &,
                                              double &,
                                              double &,
                                              double &,
                                              double &,
                                              double &);
        
        typedef void (physicsModel::*pRiemann)(const Field &,
                                               const VectorField &,
                                               const VectorField &,
                                               const VectorField &,
                                               const VectorField &,
                                               VectorField &,
                                               VectorField &,
                                               double &);
        
        typedef void (physicsModel::*pVisMod)(const Field &T,
                                              Field &miu,
                                              Field &coeff,
                                              double pt,
                                              double Mach,
                                              double gamma);
        
        int neqs;
        std::string myEoSName;
        std::string myRiemannName;
        std::string myModelName;
        std::string myVisModelName;
        
        std::map<std::string, pEoSFWD> myEosFWDMap;
        std::map<std::string, pEoSBWD> myEosBWDMap;
        std::map<std::string, pRiemann> myRiemannMap;
        std::map<std::string, pVisMod> myVisModelMap;
        std::map<std::string, std::string> myModelMap;
        
        std::map<std::string, double> myModelParameters; //parameters for EoS
    public:
        physicsModel();
        ~physicsModel();
        
        std::string getMyModelName(){return myModelName;};
        
        void solve(PDEData3d &, std::shared_ptr<EulerData3D> );
        
        int getEqNum(){return neqs;};
        
        void initial(PDEData3d &,std::shared_ptr<EulerData3D> );
    private:
        void RiemannSolver(const std::string &,
                           const Field &,
                           const VectorField &,
                           const VectorField &,
                           const VectorField &,
                           const VectorField &,
                           EulerFlux&);
        
        void con2prim(const std::string &,
                      const Field &,
                      const VectorField &,
                      VectorField &,
                      VectorField &);
        
        void prim2con(const std::string &,
                      const Field &,
                      const VectorField &,
                      VectorField &);
        
        void RiemannAUSM(
                         const Field &,
                         const VectorField &,
                         const VectorField &,
                         const VectorField &,
                         const VectorField &,
                         VectorField &,
                         VectorField &,
                         double &);
        
        void RiemannLF(
                       const Field &,
                       const VectorField &,
                       const VectorField &,
                       const VectorField &,
                       const VectorField &,
                       VectorField &,
                       VectorField &,
                       double &);
        
        void EoSIdealGasFWD(const double &,
                            const double &,
                            const double &,
                            const double &,
                            const double &, //E=p/(gamma-1)+1/2*rho*(u^2+v^2+w^2)
                            double &,
                            double &,
                            double &,
                            double &);//input rho,e output p,T
        
        void EoSIdealGasBWD(
                            const double &rho,
                            const double &u,
                            const double &v,
                            const double &w,
                            const double &p,
                            double & _rho,
                            double & _rhou,
                            double & _rhoV,
                            double & _rhow,
                            double & _E);
        
        void EoSJWLFWD(const double &,
                       const double &,
                       const double &,
                       const double &,
                       const double &, //E=p/(gamma-1)+1/2*rho*(u^2+v^2+w^2)
                       double &,
                       double &,
                       double &,
                       double &);//input rho,e output p,T
        
        void EoSJWLBWD(
                       const double &rho,
                       const double &u,
                       const double &v,
                       const double &w,
                       const double &p,
                       double & _rho,
                       double & _rhou,
                       double & _rhov,
                       double & _rhow,
                       double & _E);
        
        void sutherland(const Field &T,
                        Field &miu,
                        Field &coeff,
                        double pt,
                        double Mach,
                        double gamma);
        
        void constant(const Field &T,
                        Field &miu,
                        Field &coeff,
                        double pt,
                        double Mach,
                        double gamma);
        
        std::istream& readPhysFile(std::istream& ios);
        
        double getPressureL(const double &, const double &);
        double getPressureR(const double &, const double &);
        
        double getMachL(const double &);
        double getMachR(const double &);
        
    };
}
#endif
