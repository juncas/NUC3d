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
    BufferSend.push_back(Field(bfwidth, nz, nx));
    BufferRecv.push_back(Field(bfwidth, nz, nx));
    bufferSize[2] = mySize;
    
    BufferSend.push_back(Field(bfwidth, nz, nx));
    BufferRecv.push_back(Field(bfwidth, nz, nx));
    bufferSize[3] = mySize;
    
    mySize = nx*ny*bfwidth;
    BufferSend.push_back(Field(bfwidth, nx, ny));
    BufferRecv.push_back(Field(bfwidth, nx, ny));
    bufferSize[4] = mySize;
    
    BufferSend.push_back(Field(bfwidth, nx, ny));
    BufferRecv.push_back(Field(bfwidth, nx, ny));
    bufferSize[5] = mySize;
};

nuc3d::bufferData::~bufferData()
{
    
};
/**************************************************************************************
 Definition of member functions
 **************************************************************************************/
void nuc3d::bufferData::setBufferSend(Field& myField,int faceID)
{
    
    double *f=myField.getDataPtr();
    double *bf;
    int nx0,ny0,nz0;
    
    switch (faceID) {
        case 0:
            nx0=myField.getSizeX();
            ny0=myField.getSizeY();
            nz0=myField.getSizeZ();
            bf=BufferSend[0].getDataPtr();
            for (int k = 0; k < nz0; k++)
            {
                for (int j = 0; j < ny0; j++)
                {
                    for (int i = 0; i < bufferWidth; i++)
                    {
                        int idx_bf=bufferWidth*ny0*k+bufferWidth*j+i;
                        int idx_f=nx0*ny0*k+nx0*j+i;
                        
                        bf[idx_bf]=f[idx_f];
                    }
                    
                }
            }
            break;
        case 1:
            nx0=myField.getSizeX();
            ny0=myField.getSizeY();
            nz0=myField.getSizeZ();
            bf=BufferSend[1].getDataPtr();
            for (int k = 0; k < nz0; k++)
            {
                for (int j = 0; j < ny0; j++)
                {
                    for (int i = 0; i < bufferWidth; i++)
                    {
                        int idx_bf=bufferWidth*ny0*k+bufferWidth*j+i;
                        int idx_f=nx0*ny0*k+nx0*j+i+nx0-bufferWidth;
                        
                        bf[idx_bf]=f[idx_f];
                    }
                }
            }
            break;
        case 2:
            nx0=myField.getSizeX();
            ny0=myField.getSizeY();
            nz0=myField.getSizeZ();

            bf=BufferSend[2].getDataPtr();
            for (int k = 0; k < nz0; k++)
            {
                for (int j = 0; j < ny0; j++)
                {
                    for (int i = 0; i < bufferWidth; i++)
                    {
                        int idx_bf=bufferWidth*ny0*k+bufferWidth*j+i;
                        int idx_f=nx0*ny0*k+nx0*j+i;
                        
                        bf[idx_bf]=f[idx_f];
                    }
                    
                }
            }
            break;
            
        case 3:
            nx0=myField.getSizeX();
            ny0=myField.getSizeY();
            nz0=myField.getSizeZ();

            bf=BufferSend[3].getDataPtr();
            for (int k = 0; k < nz0; k++)
            {
                for (int j = 0; j < ny0; j++)
                {
                    for (int i = 0; i < bufferWidth; i++)
                    {
                        int idx_bf=bufferWidth*ny0*k+bufferWidth*j+i;
                        int idx_f=nx0*ny0*k+nx0*j+i+nx0-bufferWidth;
                        
                        bf[idx_bf]=f[idx_f];
                    }
                }
            }
            break;
        case 4:
            nx0=myField.getSizeX();
            ny0=myField.getSizeY();
            nz0=myField.getSizeZ();

            bf=BufferSend[4].getDataPtr();
            for (int k = 0; k < nz0; k++)
            {
                for (int j = 0; j < ny0; j++)
                {
                    for (int i = 0; i < bufferWidth; i++)
                    {
                        int idx_bf=bufferWidth*ny0*k+bufferWidth*j+i;
                        int idx_f=nx0*ny0*k+nx0*j+i;
                        
                        bf[idx_bf]=f[idx_f];
                    }
                    
                }
            }
            
            break;
        case 5:
            nx0=myField.getSizeX();
            ny0=myField.getSizeY();
            nz0=myField.getSizeZ();

            bf=BufferSend[5].getDataPtr();
            for (int k = 0; k < nz0; k++)
            {
                for (int j = 0; j < ny0; j++)
                {
                    for (int i = 0; i < bufferWidth; i++)
                    {
                        int idx_bf=bufferWidth*ny0*k+bufferWidth*j+i;
                        int idx_f=nx0*ny0*k+nx0*j+i+nx0-bufferWidth;
                        
                        bf[idx_bf]=f[idx_f];
                    }
                }
            }
            
            break;
    }
    
    
    
};
/**************************************************************************************
 End of definition
 **************************************************************************************/

