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
    protected:
        ~interoplation()=default;

    private:
        virtual void interpolationInner(
                const Field &,
                const int,
                const int,
                const int,
                const int,
                Field &) = 0;
        virtual void interpolationBoundary(
                const Field &,
                const Field &,
                const Field &,
                const int,
                const int,
                const int,
                const int,
                Field &) = 0;
    };
    
    class integration
    {
        friend class fieldOperator3d;
    protected:
        
        int nstep;
        ~integration()=default;

    private:
        virtual void getRHS(
                const Field &, // dfdx
                const Field &, // dgdy
                const Field &, // dhdz
                const double, // dt
                Field &) = 0;// Right-hand-side

        virtual void integrationAll(
                const Field &, // Right-hand-side: l*dt
                const Field &, // u(nstep)
                const Field &, // u_n
                int , // n th step
                Field &) = 0; // output
    };

    class differential
    {
        friend class fieldOperator3d;
    protected:
        ~differential()=default;

    private:
        virtual void differentialInner(
                const Field &,
                const int,
                const int,
                const int,
                Field &) = 0;

        virtual void differentialBoundary(
                const Field &,
                const Field &,
                const Field &,
                const int,
                const int,
                const int,
                Field &) = 0;
    };

    class fieldOperator3d//should only operate on fields
    {
        int bufferSize;
        
        std::vector<std::shared_ptr<interoplation>> myInteroplators;
        std::vector<std::shared_ptr<differential>> myDifferenters;

        std::shared_ptr<integration> myIntegrators;


        std::map<std::string,std::string> MethodMap;

    	public:
    		
        fieldOperator3d();
        ~fieldOperator3d();

        int getBufferSize(){return bufferSize;};
        void reconstructionInner(const Field&, int, int, Field &);

        void reconstructionBoundary(
                const Field &,
                const Field &,
                const Field &,
                const int,
                const int,
                Field &);

        void differenceInner(const Field&, int, Field &);

        void differenceBoundary(
                const Field &,
                const Field &,
                const Field &,
                const int,
                Field &);


        void timeIntegral ( const Field&, //dfdx
                            const Field&, //dfdy
                            const Field&, //dfdz
                            const Field&, // u0
                            const Field&, // un 
                            const Field&, // rhs
                            const double,
                            Field &,
                            int);
    private:
		std::istream& readIOFile(
			std::istream&,
			std::map<std::string, std::string> &);

        void setMethodIvsX();
        void setMethodIvsY();
        void setMethodIvsZ();

        void setTimeMethod(const VectorField &);

        void setDiffMethodX();
        void setDiffMethodY();
        void setDiffMethodZ();

    };
}

#endif
