//
//  bufferData.hpp
//  NUC3d
//
//  Created by Jun Peng on 15/11/5.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#ifndef bufferData_hpp
#define bufferData_hpp

#include <stdio.h>
#include "field.h"
#include "mpi.h"

namespace nuc3d
{
    class bufferData// resposible for transportation only
    {
    public:
        int bufferWidth;
        int bufferSize[6];
        
        int nx;
        int ny;
        int nz;
        
        VectorField BufferSend;
        VectorField BufferRecv;
        
        MPI_Request myRequestSend[6];
        MPI_Status myStatusSend[6];
        MPI_Request myRequestRecv[6];
        MPI_Status myStatusRecv[6];
    public:
        bufferData(int il, int jl, int kl, int bfwidth);
        ~bufferData();
        
    public:
        void setBufferSend(Field &,int);
    };
    
    typedef std::vector<bufferData> VectorBuffer;
}


#endif /* bufferData_hpp */
