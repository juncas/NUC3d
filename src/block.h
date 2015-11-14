//
//  block.hpp
//  NUC3d
//
//  Created by Jun Peng on 15/10/20.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//
#ifndef block_hpp
#define block_hpp
#include <cstdlib>
#include <memory>
#include "field.h"
#include "PDEData3d.hpp"
#include "bufferData.hpp"

namespace nuc3d
{
    class PDEData3d;
    class EulerData3D;
    class bufferData;
    class physicsModel;
    class MPIComunicator3d_nonblocking;
    class boundaryCondition;
    class IOController;
    
    
    static double lag_coeff[4][5]={
        {35.0/128.0,35.0/32.0,-35.0/64.0,7.0/32.0,-5.0/128.0},
        {-5.0/128.0,15.0/32.0,45.0/64.0,-5.0/32.0,3.0/128.0},
        {3.0/128.0,-5.0/32.0,45.0/64.0,15.0/32.0,-5.0/128.0},
        {-5.0/128.0,7.0/32.0,-35.0/64.0,35.0/32.0,35.0/128.0}
    };
    
    static double derlag_coeff[4][5]={
        {-11.0/12.0,17.0/24.0,3.0/8.0,-5.0/24.0,1.0/24.0},
        {1.0/24.0,-9.0/8.0,9.0/8.0,-1.0/24.0,0.0},
        {0.0,1.0/24.0,-9.0/8.0,9.0/8.0,-1.0/24.0},
        {-1.0/24.0,5.0/24.0,-3.0/8.0,-17.0/24.0,11.0/12.0}
    };
    
    
    class block
    {
    protected:
        int nx;
        int ny;
        int nz;
        int bfsize;
        VectorField xyz;
        VectorField xyz_center;
        
        PDEData3d myPDE;
        
        /*
         this shared point could be:
         - EulerData3D;
         - EulerReactiveData3D;
         - NaiverStokesData3d;
         - NaiverStokesReactiveData3d;
         */
        std::shared_ptr<EulerData3D> myFluxes;
        
        VectorBuffer mybuffer;
        
        double time;
        double dt;
        int istep;
        double RES;
        
    public:
        block();
        ~block();
        void initial(fieldOperator3d &,
                     physicsModel &,
                     MPIComunicator3d_nonblocking &,
                     boundaryCondition &,
                     IOController &);
        
        void solve(fieldOperator3d &,
                   physicsModel &,
                   MPIComunicator3d_nonblocking &,
                   boundaryCondition &,
                   IOController &);
        
        void printStatus();
    private:
        
       typedef double (block::*pDer_lag)(const int ,
                            const int ,
                            const int ,
                            const int ,
                            const int ,
                            const int ,
                            const Field &,
                            int,
                            int,
                            int);
        
        pDer_lag der_lag[3]={
        &nuc3d::block::interpolation_derlag_center_xi,
        &nuc3d::block::interpolation_derlag_center_eta,
        &nuc3d::block::interpolation_derlag_center_zeta};
                           
        void initialData(int,int,int,physicsModel &);
        
        void getXYZ_center();
        
        void getJacobians();
        
        void interpolation_lag(const Field &,Field &);
        
        double interpolation_lag_center(const int ,
                                        const int ,
                                        const int ,
                                        const int ,
                                        const int ,
                                        const int ,
                                        const Field &,
                                        int,
                                        int,
                                        int);
        
        void interpolation_derlag(const Field &,Field &,int);
        
        double interpolation_derlag_center_xi(const int ,
                                              const int ,
                                              const int ,
                                              const int ,
                                              const int ,
                                              const int ,
                                              const Field &,
                                              int,
                                              int,
                                              int);
        
        double interpolation_derlag_center_eta(const int ,
                                              const int ,
                                              const int ,
                                              const int ,
                                              const int ,
                                              const int ,
                                              const Field &,
                                              int,
                                              int,
                                              int);
        
        double interpolation_derlag_center_zeta(const int ,
                                              const int ,
                                              const int ,
                                              const int ,
                                              const int ,
                                              const int ,
                                              const Field &,
                                              int,
                                              int,
                                              int);
        
        void readField(std::ifstream &, Field &);
        void writeField(std::ofstream &myFile, nuc3d::Field &myField);
    };
    
}
#endif /* block_hpp */
