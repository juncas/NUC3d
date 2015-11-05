//
//  bufferData.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/11/5.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#include "bufferData.hpp"
/**************************************************************************************
 Definition of class MPIComunicatorDataLayer: contains data
 **************************************************************************************/
/**************************************************************************************
 Definition of constructors and destructors
 **************************************************************************************/
nuc3d::bufferData::bufferData(int il, int jl, int kl, int bfwidth) :
nx(il), ny(jl), nz(kl), bufferWidth(bfwidth)
{
    int mySize;
    
    mySize = bfwidth*ny*nz;
    BufferSend.push_back(Field(bfwidth, ny, nz));
    BufferRecv.push_back(Field(bfwidth, ny, nz));
    bufferSize[0] = mySize;
    
    BufferSend.push_back(Field(bfwidth, ny, nz));
    BufferRecv.push_back(Field(bfwidth, ny, nz));
    bufferSize[1] = mySize;
    
    mySize = nx*bfwidth*nz;
    BufferSend.push_back(Field(nx, bfwidth, nz));
    BufferRecv.push_back(Field(nx, bfwidth, nz));
    bufferSize[2] = mySize;
    
    BufferSend.push_back(Field(nx, bfwidth, nz));
    BufferRecv.push_back(Field(nx, bfwidth, nz));
    bufferSize[3] = mySize;
    
    mySize = nx*ny*bfwidth;
    BufferSend.push_back(Field(nx, ny, bfwidth));
    BufferRecv.push_back(Field(nx, ny, bfwidth));
    bufferSize[4] = mySize;
    
    BufferSend.push_back(Field(nx, ny, bfwidth));
    BufferRecv.push_back(Field(nx, ny, bfwidth));
    bufferSize[5] = mySize;
};

nuc3d::bufferData::~bufferData()
{
    
};
/**************************************************************************************
 Definition of member functions
 **************************************************************************************/
void nuc3d::bufferData::setBufferSend(Field& myField)
{
    for (int k = 0; k < nz; k++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < bufferWidth; i++)
            {
                BufferSend[0].setValue(i,j,k,myField.getValue(i, j, k));
                BufferSend[1].setValue(i, j, k, myField.getValue(i + (nx - bufferWidth), j, k));
            }
            
        }
    }
    
    for (int k = 0; k < nz; k++)
    {
        for (int j = 0; j < bufferWidth; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                BufferSend[2].setValue(i,j,k, myField.getValue(i, j, k));
                BufferSend[3].setValue(i, j, k, myField.getValue(i, j + (ny - bufferWidth), k));
            }
            
        }
    }
    
    for (int k = 0; k < bufferWidth; k++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                BufferSend[4].setValue(i, j, k, myField.getValue(i, j, k));
                BufferSend[5].setValue(i, j, k, myField.getValue(i, j, k + (nz - bufferWidth)));
            }
            
        }
    }
    
};
/**************************************************************************************
 End of definition
 **************************************************************************************/

