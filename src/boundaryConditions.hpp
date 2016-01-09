//
//  boundaryConditions.hpp
//  NUC3d
//
//  Created by Jun Peng on 15/11/5.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#ifndef boundaryConditions_hpp
#define boundaryConditions_hpp
#include <memory>
#include <map>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include "field.h"
#include "mpi.h"
#include "bufferData.hpp"

namespace nuc3d
{
    class EulerData3D;
    class NaiverStokesData3d;
    class MPIComunicator3d_nonblocking;
    class PDEData3d;
    class physicsModel;
    class EulerFlux;
    
    static std::string BCTypes[2]={"Exterior","Interior"};
    
    static std::string BCnames[5]={"Inlet condition",
        "Outlet condition",
        "Wall condition",
        "Symmetric condition"};
    
    class faceBC
    {
    public:
        int Type; // -1:= bc  n:= neighbour block id
        int id; // bc id or neighbour face id
        
        faceBC(int myType,int myID):Type(myType),id(myID){};
        ~faceBC(){};
    };
    
    class boundaryCondition
    {
        
        typedef void (boundaryCondition::*psetBuffer)(PDEData3d &myPDE,
                                                      physicsModel &myPhyMod,
                                                      EulerData3D &myFluxes,
                                                      VectorBuffer &myBf,
                                                      int iface);
        
        typedef void (boundaryCondition::*psetVisBuffer)(PDEData3d &myPDE,
                                                         physicsModel &myPhyMod,
                                                         NaiverStokesData3d &myFluxes,
                                                         VectorBuffer &myBf,
                                                         int iface);
        
        typedef void (boundaryCondition::*psetVisFluxBuffer)(PDEData3d &myPDE,
                                                         physicsModel &myPhyMod,
                                                         NaiverStokesData3d &myFluxes,
                                                         VectorBuffer &myBf,
                                                         int iface);
        
        psetBuffer mySetter[4]={
            &boundaryCondition::setBC_Inlet, // bc_id=0
            &boundaryCondition::setBC_Outlet,// bc_id=1
            &boundaryCondition::setBC_wall,// bc_id=2
            &boundaryCondition::setBC_symm// bc_id=3
        };
        
//        psetVisBuffer myVisSetter[4]={
//            &boundaryCondition::setVisBC_Inlet, // bc_id=0
//            &boundaryCondition::setVisBC_Outlet,// bc_id=1
//            &boundaryCondition::setVisBC_wall,// bc_id=2
//            &boundaryCondition::setVisBC_symm// bc_id=3
//        };
//
//        psetVisFluxBuffer myVisFluxSetter[4]={
//            &boundaryCondition::setVisFluxBC_Inlet, // bc_id=0
//            &boundaryCondition::setVisFluxBC_Outlet,// bc_id=1
//            &boundaryCondition::setVisFluxBC_wall,// bc_id=2
//            &boundaryCondition::setVisFluxBC_symm// bc_id=3
//        };
        
        std::vector<faceBC> BCTopo;
        
        std::vector<std::vector<double>> BCvalue;//includes necessary values for a specific BC
        
    public:
        boundaryCondition();
        ~boundaryCondition();
        
    public:
        void setBC(PDEData3d &,
                   physicsModel &,
                   EulerData3D &,
                   VectorBuffer &);
        
        void setVisBC(PDEData3d &,
                      physicsModel &,
                      NaiverStokesData3d &,
                      VectorBuffer &);
        
        void setVisFluxBC(PDEData3d &,
                      physicsModel &,
                      NaiverStokesData3d &,
                      VectorBuffer &);
        
        void initialBC(VectorBuffer &,
                       MPIComunicator3d_nonblocking &);
        
        int &getBCtype(int i){return BCTopo[i].Type;};
    private:
        std::ifstream& readBCTopo(std::ifstream&,int);
        
    //Inviscid Boundary conditions
        //intlet
        void setBC_Inlet(PDEData3d &myPDE,
                         physicsModel &myPhyMod,
                         EulerData3D &myFluxes,
                         VectorBuffer &myBf,
                         int iface);
        
        void BCsetter_inlet_xi(PDEData3d &myPDE,
                               physicsModel &myPhyMod,
                               EulerData3D &myFluxes,
                               VectorBuffer &myBf,
                               int lr);
        void BCsetter_inlet_eta(PDEData3d &myPDE,
                                physicsModel &myPhyMod,
                                EulerData3D &myFluxes,
                                VectorBuffer &myBf,
                                int lr);
        void BCsetter_inlet_zeta(PDEData3d &myPDE,
                                 physicsModel &myPhyMod,
                                 EulerData3D &myFluxes,
                                 VectorBuffer &myBf,
                                 int lr);
        
        
        //outlet
        void setBC_Outlet(PDEData3d &myPDE,
                          physicsModel &myPhyMod,
                          EulerData3D &myFluxes,
                          VectorBuffer &myBf,
                          int iface);
        
        void BCsetter_outlet_xi(PDEData3d &myPDE,
                                physicsModel &myPhyMod,
                                EulerData3D &myFluxes,
                                VectorBuffer &myBf,
                                int lr);
        void BCsetter_outlet_eta(PDEData3d &myPDE,
                                 physicsModel &myPhyMod,
                                 EulerData3D &myFluxes,
                                 VectorBuffer &myBf,
                                 int lr);
        void BCsetter_outlet_zeta(PDEData3d &myPDE,
                                  physicsModel &myPhyMod,
                                  EulerData3D &myFluxes,
                                  VectorBuffer &myBf,
                                  int lr);
        
        //wall
        void setBC_wall(PDEData3d &myPDE,
                        physicsModel &myPhyMod,
                        EulerData3D &myFluxes,
                        VectorBuffer &myBf,
                        int iface);
        
        void BCsetter_wall_xi(PDEData3d &myPDE,
                              physicsModel &myPhyMod,
                              EulerData3D &myFluxes,
                              VectorBuffer &myBf,
                              int lr);
        void BCsetter_wall_eta(PDEData3d &myPDE,
                               physicsModel &myPhyMod,
                               EulerData3D &myFluxes,
                               VectorBuffer &myBf,
                               int lr);
        void BCsetter_wall_zeta(PDEData3d &myPDE,
                                physicsModel &myPhyMod,
                                EulerData3D &myFluxes,
                                VectorBuffer &myBf,
                                int lr);
        
        //symmetric
        void setBC_symm(PDEData3d &myPDE,
                        physicsModel &myPhyMod,
                        EulerData3D &myFluxes,
                        VectorBuffer &myBf,
                        int iface);
        
        void BCsetter_symm_xi(PDEData3d &myPDE,
                              physicsModel &myPhyMod,
                              EulerData3D &myFluxes,
                              VectorBuffer &myBf,
                              int lr);
        void BCsetter_symm_eta(PDEData3d &myPDE,
                               physicsModel &myPhyMod,
                               EulerData3D &myFluxes,
                               VectorBuffer &myBf,
                               int lr);
        void BCsetter_symm_zeta(PDEData3d &myPDE,
                                physicsModel &myPhyMod,
                                EulerData3D &myFluxes,
                                VectorBuffer &myBf,
                               int lr);
   /*
        //Viscous boundary conditions
        //intlet
        void setVisBC_Inlet(PDEData3d &myPDE,
                         physicsModel &myPhyMod,
                         NaiverStokesData3d &myFluxes,
                         VectorBuffer &myBf,
                         int iface);
        
        void VisBCsetter_inlet_xi(PDEData3d &myPDE,
                               physicsModel &myPhyMod,
                               NaiverStokesData3d &myFluxes,
                               VectorBuffer &myBf,
                               int lr);
        void VisBCsetter_inlet_eta(PDEData3d &myPDE,
                                physicsModel &myPhyMod,
                                NaiverStokesData3d &myFluxes,
                                VectorBuffer &myBf,
                                int lr);
        void VisBCsetter_inlet_zeta(PDEData3d &myPDE,
                                 physicsModel &myPhyMod,
                                 NaiverStokesData3d &myFluxes,
                                 VectorBuffer &myBf,
                                 int lr);
        
        
        //outlet
        void setVisBC_Outlet(PDEData3d &myPDE,
                          physicsModel &myPhyMod,
                          NaiverStokesData3d &myFluxes,
                          VectorBuffer &myBf,
                          int iface);
        
        void VisBCsetter_outlet_xi(PDEData3d &myPDE,
                                physicsModel &myPhyMod,
                                NaiverStokesData3d &myFluxes,
                                VectorBuffer &myBf,
                                int lr);
        void VisBCsetter_outlet_eta(PDEData3d &myPDE,
                                 physicsModel &myPhyMod,
                                 NaiverStokesData3d &myFluxes,
                                 VectorBuffer &myBf,
                                 int lr);
        void VisBCsetter_outlet_zeta(PDEData3d &myPDE,
                                  physicsModel &myPhyMod,
                                  NaiverStokesData3d &myFluxes,
                                  VectorBuffer &myBf,
                                  int lr);
        
        //wall
        void setVisBC_wall(PDEData3d &myPDE,
                        physicsModel &myPhyMod,
                        NaiverStokesData3d &myFluxes,
                        VectorBuffer &myBf,
                        int iface);
        
        void VisBCsetter_wall_xi(PDEData3d &myPDE,
                              physicsModel &myPhyMod,
                              NaiverStokesData3d &myFluxes,
                              VectorBuffer &myBf,
                              int lr);
        void VisBCsetter_wall_eta(PDEData3d &myPDE,
                               physicsModel &myPhyMod,
                               NaiverStokesData3d &myFluxes,
                               VectorBuffer &myBf,
                               int lr);
        void VisBCsetter_wall_zeta(PDEData3d &myPDE,
                                physicsModel &myPhyMod,
                                NaiverStokesData3d &myFluxes,
                                VectorBuffer &myBf,
                                int lr);
        
        //symmetric
        void setVisBC_symm(PDEData3d &myPDE,
                        physicsModel &myPhyMod,
                        NaiverStokesData3d &myFluxes,
                        VectorBuffer &myBf,
                        int iface);
        
        void VisBCsetter_symm_xi(PDEData3d &myPDE,
                              physicsModel &myPhyMod,
                              NaiverStokesData3d &myFluxes,
                              VectorBuffer &myBf,
                              int lr);
        void VisBCsetter_symm_eta(PDEData3d &myPDE,
                               physicsModel &myPhyMod,
                               NaiverStokesData3d &myFluxes,
                               VectorBuffer &myBf,
                               int lr);
        void VisBCsetter_symm_zeta(PDEData3d &myPDE,
                                physicsModel &myPhyMod,
                                NaiverStokesData3d &myFluxes,
                                VectorBuffer &myBf,
                                int lr);
    //Viscous Flux boundary conditions
        //intlet
        void setVisFluxBC_Inlet(PDEData3d &myPDE,
                            physicsModel &myPhyMod,
                            NaiverStokesData3d &myFluxes,
                            VectorBuffer &myBf,
                            int iface);
        
        void VisFluxBCsetter_inlet_xi(PDEData3d &myPDE,
                                  physicsModel &myPhyMod,
                                  NaiverStokesData3d &myFluxes,
                                  VectorBuffer &myBf,
                                  int lr);
        void VisFluxBCsetter_inlet_eta(PDEData3d &myPDE,
                                   physicsModel &myPhyMod,
                                   NaiverStokesData3d &myFluxes,
                                   VectorBuffer &myBf,
                                   int lr);
        void VisFluxBCsetter_inlet_zeta(PDEData3d &myPDE,
                                    physicsModel &myPhyMod,
                                    NaiverStokesData3d &myFluxes,
                                    VectorBuffer &myBf,
                                    int lr);
        
        
        //outlet
        void setVisFluxBC_Outlet(PDEData3d &myPDE,
                             physicsModel &myPhyMod,
                             NaiverStokesData3d &myFluxes,
                             VectorBuffer &myBf,
                             int iface);
        
        void VisFluxBCsetter_outlet_xi(PDEData3d &myPDE,
                                   physicsModel &myPhyMod,
                                   NaiverStokesData3d &myFluxes,
                                   VectorBuffer &myBf,
                                   int lr);
        void VisFluxBCsetter_outlet_eta(PDEData3d &myPDE,
                                    physicsModel &myPhyMod,
                                    NaiverStokesData3d &myFluxes,
                                    VectorBuffer &myBf,
                                    int lr);
        void VisFluxBCsetter_outlet_zeta(PDEData3d &myPDE,
                                     physicsModel &myPhyMod,
                                     NaiverStokesData3d &myFluxes,
                                     VectorBuffer &myBf,
                                     int lr);
        
        //wall
        void setVisFluxBC_wall(PDEData3d &myPDE,
                           physicsModel &myPhyMod,
                           NaiverStokesData3d &myFluxes,
                           VectorBuffer &myBf,
                           int iface);
        
        void VisFluxBCsetter_wall_xi(PDEData3d &myPDE,
                                 physicsModel &myPhyMod,
                                 NaiverStokesData3d &myFluxes,
                                 VectorBuffer &myBf,
                                 int lr);
        void VisFluxBCsetter_wall_eta(PDEData3d &myPDE,
                                  physicsModel &myPhyMod,
                                  NaiverStokesData3d &myFluxes,
                                  VectorBuffer &myBf,
                                  int lr);
        void VisFluxBCsetter_wall_zeta(PDEData3d &myPDE,
                                   physicsModel &myPhyMod,
                                   NaiverStokesData3d &myFluxes,
                                   VectorBuffer &myBf,
                                   int lr);
        
        //symmetric
        void setVisFluxBC_symm(PDEData3d &myPDE,
                           physicsModel &myPhyMod,
                           NaiverStokesData3d &myFluxes,
                           VectorBuffer &myBf,
                           int iface);
        
        void VisFluxBCsetter_symm_xi(PDEData3d &myPDE,
                                 physicsModel &myPhyMod,
                                 NaiverStokesData3d &myFluxes,
                                 VectorBuffer &myBf,
                                 int lr);
        void VisFluxBCsetter_symm_eta(PDEData3d &myPDE,
                                  physicsModel &myPhyMod,
                                  NaiverStokesData3d &myFluxes,
                                  VectorBuffer &myBf,
                                  int lr);
        void VisFluxBCsetter_symm_zeta(PDEData3d &myPDE,
                                   physicsModel &myPhyMod,
                                   NaiverStokesData3d &myFluxes,
                                   VectorBuffer &myBf,
                                   int lr);
    */

    };
}

#endif /* boundaryConditions_hpp */
