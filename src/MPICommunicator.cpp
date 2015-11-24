#ifndef MPICommunicator3d_cpp
#define MPICommunicator3d_cpp
#include "MPICommunicator.h"
#include "bufferData.hpp"
#define BoundaryFace (-1)
/**************************************************************************************
	Definition of class MPIComunicator3d_nonblocking: non_blocking coommunicator
 **************************************************************************************/
/**************************************************************************************
 Definition of constructors and destructors
 **************************************************************************************/
nuc3d::MPIComunicator3d_nonblocking::MPIComunicator3d_nonblocking() :
MPICommunicator()
{
    NeighBourFaces = new int[6];
    NeighBours = new int[6];
};

nuc3d::MPIComunicator3d_nonblocking::~MPIComunicator3d_nonblocking()
{
    delete [] NeighBours;
    delete[] NeighBourFaces;
};
/**************************************************************************************
 Definition of member functions
 **************************************************************************************/
void nuc3d::MPIComunicator3d_nonblocking::setTopo(std::vector<faceBC> &bcTopo)
{
    for (auto iter = bcTopo.begin(); iter !=bcTopo.end(); iter++)
    {
        NeighBours[iter-bcTopo.begin()] = ((iter->Type != BoundaryFace) ? iter->Type : MPI_PROC_NULL);
        NeighBourFaces[iter-bcTopo.begin()] = iter->id;
    }
};

void nuc3d::MPIComunicator3d_nonblocking::bufferSendRecv(Field &myfield, bufferData& buffer,int dir,int fdID)
{
    switch(dir)
    {
        case 0:
            bufferSendRecvFace(myfield,buffer,0,fdID);
            bufferSendRecvFace(myfield,buffer,1,fdID);
            break;
        case 1:
            bufferSendRecvFace(myfield,buffer,2,fdID);
            bufferSendRecvFace(myfield,buffer,3,fdID);
            break;
        case 2:
            bufferSendRecvFace(myfield,buffer,4,fdID);
            bufferSendRecvFace(myfield,buffer,5,fdID);
            break;
    }
};

void nuc3d::MPIComunicator3d_nonblocking::bufferSendRecvFace(Field &myfield, bufferData &buffer,int i,int fdID)
{
    if(NeighBours[i]!=MPI_PROC_NULL)
    {
        buffer.setBufferSend(myfield, i);
        int sendTag = setMPISendTag(this->getMyId(), i ,fdID);
        int RecvTag = setMPISendTag(NeighBours[i], NeighBourFaces[i],fdID);
        
        MPI_Irecv(buffer.BufferRecv[i].getDataPtr(),
                  buffer.bufferSize[i],
                  MPI_DOUBLE,
                  NeighBours[i],
                  RecvTag,
                  MPI_COMM_WORLD,
                  &(buffer.myRequestRecv[i]));
        
        MPI_Isend(buffer.BufferSend[i].getDataPtr(),
                  buffer.bufferSize[i],
                  MPI_DOUBLE,
                  NeighBours[i],
                  sendTag,
                  MPI_COMM_WORLD,
                  &(buffer.myRequestSend[i]));
    }
    
    
}

void nuc3d::MPIComunicator3d_nonblocking::waitSendRecv(bufferData& buffer,int dir)
{
    switch(dir)
    {
        case 0:
            waitSendRecvFace(buffer,0);
            waitSendRecvFace(buffer,1);
            break;
        case 1:
            waitSendRecvFace(buffer,2);
            waitSendRecvFace(buffer,3);
            break;
        case 2:
            waitSendRecvFace(buffer,4);
            waitSendRecvFace(buffer,5);
            break;
    }
    
};
void nuc3d::MPIComunicator3d_nonblocking::waitSendRecvFace(bufferData& buffer,int bfID)
{
    if(NeighBours[bfID]!=MPI_PROC_NULL)
    {
        if((MPI_Wait(&buffer.myRequestSend[bfID], &buffer.myStatusSend[bfID])!= MPI_SUCCESS)
           ||(MPI_Wait(&buffer.myRequestRecv[bfID], &buffer.myStatusRecv[bfID])!= MPI_SUCCESS))
        {
            std::cout << "Communicating failed" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 999);
        }
    }
    
    
}


/**************************************************************************************
 End of definition
 **************************************************************************************/
#endif
