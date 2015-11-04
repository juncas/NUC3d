//
//  PDEData3d.hpp
//  NUC3d
//
//  Created by Jun Peng on 15/11/4.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#ifndef PDEData3d_hpp
#define PDEData3d_hpp
#include <stdio.h>
#include "fieldOperator.h"

namespace nuc3d
{
    class EulerData3D;
    class EulerFlux;
    class PDEData3d;
    class physicsModel;
    
    class PDEData3d
    {
        int nEquations;
        
        //a PDE is consist of such four vectors
        // other vectors are intermediate data
        /*
         
         dQ
         -- =   -RHS
         dt
         
         
         df     dg     dh
         RHS = ---- + ---- + ----- - source
         dxi    deta   dzeta
         */
        
        //current vector
        VectorField Q_Euler;
        VectorField Q_work;
        
        VectorField RHS;
        
        double dt_local;
        double dt_global;
        double res_local;
        double res_global;
        
    public:
        PDEData3d();
        PDEData3d(int,int,int,int);
        ~PDEData3d();
        
        
        void initPDEData3d(int,int,int,int);
        
        VectorField& getRHS();
        VectorField& getQ();
        void setDt(double);
        
        void solve(fieldOperator3d &,
                   int);
        
        void setRES();
    };
}

#endif /* PDEData3d_hpp */
