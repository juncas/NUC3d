#ifndef physicsModel_h
#define physicsModel_h

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <cmath>
#include <memory>
#include "field.h"

namespace nuc3d
{
    class EulerData3D;
    class EulerFlux;
    class PDEData3d;
    
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
        
        //        typedef void (physicsModel::*pRiemann)(const Field &,
        //                                               const VectorField &,
        //                                               const VectorField &,
        //                                               const VectorField &,
        //                                               const VectorField &,
        //                                               VectorField &,
        //                                               VectorField &,
        //                                               double &);
        
        typedef void (physicsModel::*pRiemann)(const double &U0,
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
                                               double *fluxn);
        
        typedef void (physicsModel::*pRiemannPoint)(const std::vector<double> &,
                                                    const double ,
                                                    const double ,
                                                    const double ,
                                                    const double ,
                                                    std::vector<double> &,
                                                    std::vector<double> &);
        
        typedef void (physicsModel::*pVisMod)(const Field &T,
                                              Field &miu,
                                              Field &coeff,
                                              double Reynolds,
                                              double Mach,
                                              double pt,
                                              double gamma,
                                              double T_ref,
                                              double T_inf);
        
        int neqs;
        std::string myEoSName;
        std::string myRiemannName;
        std::string myModelName;
        std::string myVisModelName;
        
        std::map<std::string, pEoSFWD> myEosFWDMap;
        std::map<std::string, pEoSBWD> myEosBWDMap;
        std::map<std::string, pRiemann> myRiemannMap;
        std::map<std::string, pVisMod> myVisModelMap;
        std::map<std::string, pRiemannPoint> myRiemannPointMap;
        
        std::map<std::string, double> myModelParameters; //parameters for EoS
    public:
        physicsModel();
        ~physicsModel();
        
        std::string getMyModelName(){return myModelName;};
        
        void solve(PDEData3d &, EulerData3D *);
        
        void solveRiemann(PDEData3d &myPDE,EulerData3D *myEuler);
        
        void solveRiemannPoint(const std::vector<double> &prim,
                               const double jac,
                               const double xx_x,
                               const double xx_y,
                               const double xx_z,
                               std::vector<double> &fluxl,
                               std::vector<double> &fluxr);
        
        
        void getMiu(Field &,
                    Field &,
                    Field &);
        
        int getEqNum(){return neqs;};
        double getMach(){return myModelParameters["Mach"];};
        double getGamma(){return myModelParameters["Gamma"];};
        double getWallTemp(){return myModelParameters["T_wall"];};
        
        void initial(PDEData3d &,std::shared_ptr<EulerData3D> );
        
        void prim2conPoint(const double &rho,
                           const double &u,
                           const double &v,
                           const double &w,
                           const double &p,
                           double &T,
                           double &e,
                           double &alpha);
        
        void getPrim(Field &jac,VectorField &Q,VectorField &Prim, VectorField &Acoust);
    private:
//        void RiemannSolver(const std::string &,
//                           const Field &,
//                           const VectorField &,
//                           const VectorField &,
//                           const VectorField &,
//                           const VectorField &,
//                           EulerFlux&);
        
        void RiemannSolver(const std::string &SolverName,
                           const Field &Jacobian,
                           const VectorField &xi_xyz,
                           const VectorField &eta_xyz,
                           const VectorField &zeta_xyz,
                           const VectorField &W_vec,
                           const VectorField &W0_vec,
                           const VectorField &Q_vec,
                           EulerFlux& myFlux_xi,
                           EulerFlux& myFlux_eta,
                           EulerFlux& myFlux_zeta);
        
        void con2prim(const std::string &,
                      const Field &,
                      VectorField &,
                      VectorField &,
                      VectorField &);
        
        void prim2con(const std::string &,
                      const Field &,
                      const VectorField &,
                      VectorField &);
        
        void RiemannAUSMp(const double &U0,
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
                         double *fluxn);
        void RiemannAUSM(const double &U0,
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
                         double *fluxn);
        
        void RiemannLF(const double &U0,
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
                         double *fluxn);
        
        void EoSIdealGasFWD(const double &,
                            const double &,
                            const double &,
                            const double &,
                            const double &, //E=p/(gamma-1)+1/2*rho*(u^2+v^2+w^2)
                            double &,
                            double &,
                            double &,
                            double &);//input rho,e output p,T
        
        void EoSIdealGasBWD(const double &rho,
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
        
        void EoSJWLBWD(const double &rho,
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
                        double Reynolds,
                        double Mach,
                        double pt,
                        double gamma,
                        double T_ref,
                        double T_inf);
        
        void constant(const Field &T,
                      Field &miu,
                      Field &coeff,
                      double Reynolds,
                      double Mach,
                      double pt,
                      double gamma,
                      double T_ref,
                      double T_inf);
        
        std::istream& readPhysFile(std::istream& ios);
        
        double getPressureL(const double &, const double &);
        double getPressureR(const double &, const double &);
        
        double getMachL(const double &);
        double getMachR(const double &);
        
        double getPressureLp(const double &, const double &);
        double getPressureRp(const double &, const double &);
        
        double getMachLp(const double &);
        double getMachRp(const double &);
        
        double EoSIdealGasgetPressure(const double &rho,
                                      const double &u,
                                      const double &v,
                                      const double &w,
                                      const double &E,
                                      const double &mach,
                                      const double &gamma);
        
        double EoSIdealGasgetTemp(const double &rho,
                                  const double &u,
                                  const double &v,
                                  const double &w,
                                  const double &E,
                                  const double &p,
                                  const double &mach,
                                  const double &gamma);
        
        double EoSIdealGasgetE0(const double &rho,
                                const double &u,
                                const double &v,
                                const double &w,
                                const double &E,
                                const double &p,
                                const double &mach,
                                const double &gamma);
        
        double EoSIdealGasgetAlpha(const double &rho,
                                   const double &u,
                                   const double &v,
                                   const double &w,
                                   const double &E,
                                   const double &p,
                                   const double &mach,
                                   const double &gamma);
        
        void solveRiemannPointAUSM(const std::vector<double> &prim,
                                   const double jac,
                                   const double xx_x,
                                   const double xx_y,
                                   const double xx_z,
                                   std::vector<double> &fluxl,
                                   std::vector<double> &fluxr);
        
        void solveRiemannPointAUSMp(const std::vector<double> &prim,
                                   const double jac,
                                   const double xx_x,
                                   const double xx_y,
                                   const double xx_z,
                                   std::vector<double> &fluxl,
                                   std::vector<double> &fluxr);
        
        void solveRiemannPointLF(const std::vector<double> &prim,
                                 const double jac,
                                 const double xx_x,
                                 const double xx_y,
                                 const double xx_z,
                                 std::vector<double> &fluxl,
                                 std::vector<double> &fluxr);
        
        
    };
}
#endif
