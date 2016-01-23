#ifndef fieldOperator3d_h
#define fieldOperator3d_h

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <cmath>
#include <memory>
#include "field.h"

namespace nuc3d
{
    class interoplation
    {
        friend class fieldOperator3d;
    public:
        interoplation(){};
        ~interoplation(){};
        
    private:
        virtual void interpolationInner(const Field &,
                                        const int,
                                        Field &,
                                        const int)=0;
        virtual void interpolationBoundaryL(
                                            const Field &,
                                            const Field &,
                                            const int,
                                            Field &,
                                            const int)=0;
        virtual void interpolationBoundaryR(
                                            const Field &,
                                            const Field &,
                                            const int,
                                            Field &,
                                            const int)=0;
    };
    
    class integration
    {
        friend class fieldOperator3d;
    protected:
        
        int nstep;
        
    private:
        
        virtual void integrationAll(const VectorField &, // Right-hand-side: l*dt
                                    const VectorField &, // u_n
                                    VectorField &, // u_i
                                    double,
                                    int)=0; // step n
    public:
        integration(int i){nstep=i;};
        ~integration()=default;
    };
    
    class differential
    {
        friend class fieldOperator3d;
    public:
        inline differential(){};
        inline ~differential(){};
        
    private:
        virtual void differentialInner(const Field &,
                                       Field &,
                                       const int)=0;
        
        virtual void differentialBoundaryL(const Field &,
                                           const Field &,
                                           Field &,
                                           const int)=0;
        virtual void differentialBoundaryR(const Field &,
                                           const Field &,
                                           Field &,
                                           const int)=0;
    };
    
    class fieldOperator3d//should only operate on fields
    {
        
        std::vector<std::shared_ptr<interoplation>> myInteroplators;
        std::vector<std::shared_ptr<interoplation>> myInteroplatorsBND;
        std::vector<std::shared_ptr<differential>> myDifferenters;
        std::vector<std::shared_ptr<differential>> myDifferentersBND;
        
        std::shared_ptr<integration> myIntegrators;
        
        
        std::map<std::string,std::string> MethodMap;
        
    public:
        
        int bufferSize_xi;
        int bufferSize_eta;
        int bufferSize_zeta;
        
        fieldOperator3d();
        ~fieldOperator3d();
        void initial(const VectorField &);
        
        int getBufferSize();
        int getSteps(){return myIntegrators->nstep;}
        
        void reconstructionInner(const Field&, int, int, Field &);
        
        void reconstructionBoundary(const Field &,
                                    const Field &,
                                    const Field &,
                                    int,
                                    int,
                                    Field &,
                                    int typeL,
                                    int typeR);
        
        void differenceInner(const Field&, int, Field &);
        
        void differenceBoundary(
                                const Field &,
                                const Field &,
                                const Field &,
                                int,
                                Field &,
                                int typeL,
                                int typeR);
        
        
        void timeIntegral (const VectorField&, // rhs
                           const VectorField&, // u_n
                           VectorField &, //u_i
                           double,
                           int);
    private:
        std::istream& readIOFile(
                                 std::istream&,
                                 std::map<std::string, std::string> &);
        
        void setMethodIvsX();
        void setMethodIvsY();
        void setMethodIvsZ();
        
        void setTimeMethod();
        
        void setDiffMethodX();
        void setDiffMethodY();
        void setDiffMethodZ();
        
        void reconstructionBoundaryInnerL(const Field &fieldIN,
                                          const Field &boundaryL,
                                          const int direction,
                                          const int upwind,
                                          Field &fieldOUT);
        
        void reconstructionBoundaryInnerR(const Field &fieldIN,
                                          const Field &boundaryR,
                                          const int direction,
                                          const int upwind,
                                          Field &fieldOUT);
        
        
        void reconstructionBoundaryExteriorL(const Field &fieldIN,
                                             const Field &boundaryL,
                                             const int direction,
                                             const int upwind,
                                             Field &fieldOUT);
        
        void reconstructionBoundaryExteriorR(const Field &fieldIN,
                                             const Field &boundaryR,
                                             const int direction,
                                             const int upwind,
                                             Field &fieldOUT);
        
        void differenceBoundaryInnerL(const Field &fieldIN,
                                      const Field &boundaryL,
                                      const int direction,
                                      Field &fieldOUT);
        
        void differenceBoundaryExteriorL(const Field &fieldIN,
                                         const Field &boundaryL,
                                         const int direction,                                             Field &fieldOUT);
        void differenceBoundaryInnerR(const Field &fieldIN,
                                      const Field &boundaryR,
                                      const int direction,
                                      Field &fieldOUT);
        
        void differenceBoundaryExteriorR(const Field &fieldIN,
                                         const Field &boundaryR,
                                         const int direction,
                                         Field &fieldOUT);
        
        
    };
}

#endif
