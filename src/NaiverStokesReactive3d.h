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
    virtual void solveLocal();
    virtual void solve(fieldOperator3d &,bufferData &);
    virtual void solveInv(fieldOperator3d &,bufferData &);
    virtual void solveVis(fieldOperator3d &,bufferData &);
    virtual void solveSource();
};
}
#endif /* NaiverStokesReactive3d_hpp */
