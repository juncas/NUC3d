#ifndef MPICommunicator3d_h
#define MPICommunicator3d_h
# include "field.h"
#include "MPICommunicatorBase.h"
#include "boundaryConditions.hpp"

namespace nuc3d
{
    class bufferData;
    
    class MPIComunicator3d_nonblocking : public MPICommunicator
	{
		int* NeighBourFaces;
		int* NeighBours;
		int* sendTags;
		int* recvTags;

	public:
		MPIComunicator3d_nonblocking();
		~MPIComunicator3d_nonblocking();

	public:
		void setTopo(std::vector<faceBC> &);
		void bufferSendRecv(bufferData &,int);
		void waitAllSendRecv(bufferData& buffer);
		void ReduceData(double );
	private:
		int setMPISendTag(
			int ProcId,
			int FaceId
			)
		{
			return 1 + FaceId * 10000 + ProcId * 100;
		};
		int setMPIRecvTag(
			int ProcId,
			int FaceId
			)
		{
			return  FaceId * 10000 + ProcId * 100;
		};

	};
}
#endif