//
//  eulerReactive3d.hpp
//  NUC3d
//
//  Created by Jun Peng on 15/11/3.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#ifndef eulerReactive3d_hpp
#define eulerReactive3d_hpp
#include "euler3d.h"

namespace nuc3d
{
    class EulerReactiveData3D : virtual public EulerData3D
    {
        friend class physicsModel;
        
    protected:
        VectorField source_xi;
        VectorField source_eta;
        VectorField source_zeta;
        
    public:
        EulerReactiveData3D(int,int,int,int);
        ~EulerReactiveData3D();
        
        virtual VectorField& getDrivativeXi();
        virtual VectorField& getDrivativeEta();
        virtual VectorField& getDrivativeZeta();
        
        virtual void solve(PDEData3d &,
                           fieldOperator3d &,
                           bufferData &,
                           physicsModel &,
                           MPIComunicator3d_nonblocking &);
        
    protected:
        void solveSource(physicsModel &);
        void setSource();
    private:
        virtual void solveRHS(PDEData3d &);

    };
}

#endif /* eulerReactive3d_hpp */
