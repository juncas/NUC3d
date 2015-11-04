//
//  NaiverStokesReactive3d.hpp
//  NUC3d
//
//  Created by Jun Peng on 15/11/3.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#ifndef NaiverStokesReactive3d_hpp
#define NaiverStokesReactive3d_hpp
#include "eulerReactive3d.h"
#include "NaiverStokes3d.h"

namespace nuc3d
{
class NaiverStokesReactiveData3d : public EulerReactiveData3D, public NaiverStokesData3d
{
public:
    NaiverStokesReactiveData3d(int,int,int,int);
    ~NaiverStokesReactiveData3d();
    
    virtual VectorField& getDrivativeXi();
    virtual VectorField& getDrivativeEta();
    virtual VectorField& getDrivativeZeta();
    
    virtual void solve(PDEData3d &,
                       fieldOperator3d &,
                       bufferData &,
                       physicsModel &,
                       MPIComunicator3d_nonblocking &);
private:
    virtual void solveRHS(PDEData3d &);
};
}
#endif /* NaiverStokesReactive3d_hpp */
