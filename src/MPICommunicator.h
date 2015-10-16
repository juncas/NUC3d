#ifndef MPICommunicator3d_h
#define MPICommunicator3d_h
# include "field.h"
#include "MPICommunicatorBase.h"

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
		void setBufferSend(Field &);
		void deallocatebuffer();
	};


	class MPIComunicator3d_nonblocking : public MPICommunicator
	{
		int* NeighBours;
		int* NeighBourFaces;
		int* NeighBourBlocks;
		int* sendTags;
		int* recvTags;

	public:
		MPIComunicator3d_nonblocking(int&, char **&);
		~MPIComunicator3d_nonblocking();

	public:
		void setTopo(int *, int *,int *);
		void bufferSendRecv(bufferData &,int);
		void waitAllSendRecv(bufferData& buffer);
		void ReduceData(double );
	private:
		int setMPISendTag(
			int ProcId,
			int BlockId,
			int FaceId
			)
		{
			return 1 + FaceId * 10 + BlockId * 100 + ProcId * 10000;
		};
		int setMPIRecvTag(
			int ProcId,
			int BlockId,
			int FaceId
			)
		{
			return  FaceId * 10 + BlockId * 100 + ProcId * 10000;
		};

	};
}
#endif