#ifndef MPICommunicator3d_h
#define MPICommunicator3d_h
# include "field.h"
#include "MPICommunicatorBase.h"

namespace nuc3d
{
    class bufferData;
    
    class MPIComunicator3d_nonblocking : public MPICommunicator
	{
		int* NeighBours;
		int* NeighBourFaces;
		int* NeighBourBlocks;
		int* sendTags;
		int* recvTags;

	public:
		MPIComunicator3d_nonblocking();
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