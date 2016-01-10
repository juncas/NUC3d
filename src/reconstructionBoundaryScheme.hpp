//
//  reconstructionBoundaryScheme.hpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/10.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//

#ifndef reconstructionBoundaryScheme_hpp
#define reconstructionBoundaryScheme_hpp

#include <stdio.h>

#include "fieldOperator.h"

namespace nuc3d
{
    
    class ReconstructionboundaryScheme: public interoplation
    {
    public:
        ReconstructionboundaryScheme();
        ~ReconstructionboundaryScheme();
        // 6th order central difference functions
        
        void interpolationInner(const Field &,
                                const int,
                                Field &,
                                const int);
        
        void interpolationBoundaryL(const Field &,
                                    const Field &,
                                    const int,
                                    Field &,
                                    const int);
        
        void interpolationBoundaryR(const Field &,
                                    const Field &,
                                    const int,
                                    Field &,
                                    const int);
    private:
        void boundary_pL(const Field & fieldIN,
                         const Field & boundaryL,
                         Field & fieldOUT,
                         const int tilesize);
        
        void boundary_pR(const Field & fieldIN,
                         const Field & boundaryL,
                         Field & fieldOUT,
                         const int tilesize);
        
        void firstOrderP(double *,double *);
        void secondOrderP(double *,double *);
        void thirdOrderP(double *,double *);
        void fifthOrderP(double *,double *);
        
        void boundary_nL(const Field & fieldIN,
                         const Field & boundaryL,
                         Field & fieldOUT,
                         const int tilesize);
        void boundary_nR(const Field & fieldIN,
                         const Field & boundaryL,
                         Field & fieldOUT,
                         const int tilesize);
        
        void firstOrderN(double *,double *);
        void secondOrderN(double *,double *);
        void thirdOrderN(double *,double *);
        void fifthOrderN(double *,double *);
        
    };
}
#endif /* reconstructionBoundaryScheme_hpp */
