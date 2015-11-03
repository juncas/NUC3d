#ifndef MPICommunicator3d_cpp
#define MPICommunicator3d_cpp
#include "MPICommunicator.h"
#define BoundaryFace (-1)
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


void nuc3d::bufferData::deallocatebuffer()
{
}
/**************************************************************************************
								End of definition
**************************************************************************************/

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