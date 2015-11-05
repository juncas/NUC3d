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
	NeighBours = new int[6];
	NeighBourFaces = new int[6];
	NeighBourBlocks = new int[6];
};

nuc3d::MPIComunicator3d_nonblocking::~MPIComunicator3d_nonblocking()
{
	delete[] NeighBours;
	delete[] NeighBourFaces;
};
/**************************************************************************************
						Definition of member functions
**************************************************************************************/
void nuc3d::MPIComunicator3d_nonblocking::setTopo(int *myNeibours, int *myNeighbourFaces, int *myNeighBlocks)
{
	for (int i = 0; i < 6; i++)
	{
		NeighBours[i] = ((myNeibours[i] != BoundaryFace) ? myNeibours[i] : MPI_PROC_NULL);
		NeighBourFaces[i] = myNeighbourFaces[i];
	}
};

void nuc3d::MPIComunicator3d_nonblocking::bufferSendRecv(bufferData& buffer, int blockId)
{
	for (int i = 0; i < 6; i++)
	{
		int sendTag = setMPISendTag(this->getMyId(), blockId, i);
		int RecvTag = setMPISendTag(NeighBours[i], NeighBourBlocks[i],NeighBourFaces[i]);

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
};

void nuc3d::MPIComunicator3d_nonblocking::waitAllSendRecv(bufferData& buffer)
{
	if ((MPI_Waitall(6, buffer.myRequestSend, buffer.myStatusSend) != MPI_SUCCESS) ||
		(MPI_Waitall(6, buffer.myRequestRecv, buffer.myStatusRecv) != MPI_SUCCESS))
	{
		std::cout << "Communicating failed" << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 999);
	}
};

/**************************************************************************************
								End of definition
**************************************************************************************/
#endif