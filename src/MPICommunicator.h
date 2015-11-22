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
        void bufferSendRecv(Field &, bufferData &,int,int);
        void waitSendRecv(bufferData& buffer,int);
        void ReduceData(double );
        inline int getFaceType(int facdID){ return NeighBours[facdID];};
    private:
        void bufferSendRecvFace(Field &, bufferData &,int,int);
        void waitSendRecvFace(bufferData& buffer,int);
        
        
        int setMPISendTag(
                          int ProcId,
                          int FaceId,
                          int FieldID)
        {
            return FieldID + FaceId * 100 + ProcId * 1000;
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