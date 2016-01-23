#ifndef fieldOperator3d_cpp
#define fieldOperator3d_cpp

#include "fieldOperator.h"
#include "weno7js.hpp"
#include "weno5js.h"
#include "weno5z.hpp"
#include "weno3js.hpp"
#include "crweno5.hpp"
#include "upwind1st.hpp"
#include "TVD-RK3.h"
#include "centraldifference2nd.hpp"
#include "centraldifference6th.h"
#include "centraldifference8th.hpp"
#include "boundaryScheme.hpp"
#include "reconstructionBoundaryScheme.hpp"

nuc3d::fieldOperator3d::fieldOperator3d():
MethodMap(
          {
              {"scheme_time","RK3"},
              {"scheme_x_ivs","weno5js"},
              {"scheme_y_ivs","weno5js"},
              {"scheme_z_ivs","weno5js"},
              {"scheme_x_vis","cd6"},
              {"scheme_y_vis","cd6"},
              {"scheme_z_vis","cd6"}
          }
          )
{
    std::ifstream file("inp/Method.in");
    auto count=MethodMap.size();
    //std::cout<<"Reading Method input file..."<<std::endl;
    if(file)
    {
        while(readIOFile(file,MethodMap))
            count--;
        if(count!=0)
        {
            std::cout<<"Error, Missing IOController parameters in file IOController.io!"<<std::endl;
            std::cout<<"Using Default value!"<<std::endl;
        }
        
    }
    else
    {
        std::cout<<"IO file \'IOController.in\' does not exist!"
        <<std::endl;
        exit(0);
    }
    
    file.close();
    
    setMethodIvsX();
    setMethodIvsY();
    setMethodIvsZ();
    
    setDiffMethodX();
    setDiffMethodY();
    setDiffMethodZ();
    
    setTimeMethod();
    
    if( myInteroplators.size()!=3)
    {
        std::cout<<"Number of Reconstruction method: "<<myInteroplators.size()
        <<std::endl;
    }
    
    
    if(myDifferenters.size()!=3)
    {
        std::cout<<"Number of Difference method:"<<myDifferenters.size()
        <<std::endl;
    }
    
    
    if(myIntegrators==nullptr)
    {
        std::cout<<"Initial time method failed!"
        <<std::endl;
    }
    //std::cout<<"Method.in has been read!"<<std::endl;
    
}

nuc3d::fieldOperator3d::~fieldOperator3d()
{
    
}
/**************************************************************************************
 Definition of member functions
 **************************************************************************************/
/**************************************************************************************
 I. General IO functions
 **************************************************************************************/
std::istream& nuc3d::fieldOperator3d::readIOFile(
                                                 std::istream& ios,
                                                 std::map<std::string,std::string> &Methods)
{
    std::string word0;
    std::string word1;
    ios>>word0>>word1;
    
    
    if(Methods.find(word0)!=Methods.end())
    {
        std::istringstream value0(word1);
        std::string value;
        value0>>value;
        Methods[word0]=value;
    }
    else if(word0.size())
    {
        std::cout<<"Parameter \'"<<word0<<"in file \'IOController.io\' is not supported!"<<std::endl;
        std::cout<<"Make sure that every parameter exist in the following name list:"
        <<std::endl;
        for(auto beg=Methods.begin();beg!=Methods.end();beg++)
            std::cout<<beg->first<<" "<<beg->second<<std::endl;
    }
    return ios;
}


void nuc3d::fieldOperator3d::setMethodIvsX()
{
    std::string s=MethodMap["scheme_x_ivs"];
    
    if(s=="weno7js")
    {
        bufferSize_xi=4;
        myInteroplators.push_back(std::make_shared<weno7js>());
    }
    if(s=="weno5js")
    {
        bufferSize_xi=3;
        myInteroplators.push_back(std::make_shared<weno5js>());
    }
    else if(s=="weno5z")
    {
        bufferSize_xi=3;
        myInteroplators.push_back(std::make_shared<weno5z>());
    }
    else if(s=="crweno5")
    {
        bufferSize_xi=3;
        myInteroplators.push_back(std::make_shared<crweno5>());
    }
    else if(s=="weno3js")
    {
        bufferSize_xi=2;
        myInteroplators.push_back(std::make_shared<weno3js>());
    }
    else if(s=="upwind1st")
    {
        bufferSize_xi=1;
        myInteroplators.push_back(std::make_shared<upwind1st>());
    }
    else
    {
        myInteroplators.push_back(std::make_shared<weno5js>());
    }
    
    myInteroplatorsBND.push_back(std::make_shared<ReconstructionboundaryScheme>());
}


void nuc3d::fieldOperator3d::setMethodIvsY()
{
    std::string s=MethodMap["scheme_y_ivs"];
    
    if(s=="weno7js")
    {
        bufferSize_eta=4;
        myInteroplators.push_back(std::make_shared<weno7js>());
    }    if(s=="weno5js")
    {
        bufferSize_eta=3;
        myInteroplators.push_back(std::make_shared<weno5js>());
    }
    else if(s=="weno5z")
    {
        bufferSize_eta=3;
        myInteroplators.push_back(std::make_shared<weno5z>());
    }
    else if(s=="crweno5")
    {
        bufferSize_eta=3;
        myInteroplators.push_back(std::make_shared<crweno5>());
    }
    else if(s=="weno3js")
    {
        bufferSize_eta=2;
        myInteroplators.push_back(std::make_shared<weno3js>());
    }
    else if(s=="upwind1st")
    {
        bufferSize_eta=1;
        myInteroplators.push_back(std::make_shared<upwind1st>());
    }
    else
    {
        myInteroplators.push_back(std::make_shared<weno5js>());
    }
    myInteroplatorsBND.push_back(std::make_shared<ReconstructionboundaryScheme>());
    
}

void nuc3d::fieldOperator3d::setMethodIvsZ()
{
    std::string s=MethodMap["scheme_z_ivs"];
    
    if(s=="weno7js")
    {
        bufferSize_zeta=4;
        myInteroplators.push_back(std::make_shared<weno7js>());
    }
    if(s=="weno5js")
    {
        bufferSize_zeta=3;
        myInteroplators.push_back(std::make_shared<weno5js>());
    }
    else if(s=="weno5z")
    {
        bufferSize_zeta=3;
        myInteroplators.push_back(std::make_shared<weno5z>());
    }
    else if(s=="crweno5")
    {
        bufferSize_zeta=3;
        myInteroplators.push_back(std::make_shared<crweno5>());
    }
    else if(s=="weno3js")
    {
        bufferSize_zeta=2;
        myInteroplators.push_back(std::make_shared<weno3js>());
    }
    else if(s=="upwind1st")
    {
        bufferSize_zeta=1;
        myInteroplators.push_back(std::make_shared<upwind1st>());
    }
    else
    {
        myInteroplators.push_back(std::make_shared<weno5js>());
    }
    myInteroplatorsBND.push_back(std::make_shared<ReconstructionboundaryScheme>());
    
    
}

void nuc3d::fieldOperator3d::setDiffMethodX()
{
    std::string s=MethodMap["scheme_x_vis"];
    
    if(s=="cd8")
    {
        myDifferenters.push_back(std::make_shared<centraldifference8th>());
    }
    if(s=="cd6")
    {
        myDifferenters.push_back(std::make_shared<centraldifference6th>());
    }
    else if(s=="cd2")
    {
        myDifferenters.push_back(std::make_shared<centraldifference2nd>());
    }
    else
    {
        myDifferenters.push_back(std::make_shared<centraldifference6th>());
    }
    
    myDifferentersBND.push_back(std::make_shared<boundaryScheme>());
    
}

void nuc3d::fieldOperator3d::setDiffMethodY()
{
    std::string s=MethodMap["scheme_y_vis"];
    
    if(s=="cd8")
    {
        myDifferenters.push_back(std::make_shared<centraldifference8th>());
    }
    if(s=="cd6")
    {
        myDifferenters.push_back(std::make_shared<centraldifference6th>());
    }
    else if(s=="cd2")
    {
        myDifferenters.push_back(std::make_shared<centraldifference2nd>());
    }
    else
    {
        myDifferenters.push_back(std::make_shared<centraldifference6th>());
    }
    
    myDifferentersBND.push_back(std::make_shared<boundaryScheme>());
}

void nuc3d::fieldOperator3d::setDiffMethodZ()
{
    std::string s=MethodMap["scheme_z_vis"];
    
    if(s=="cd8")
    {
        myDifferenters.push_back(std::make_shared<centraldifference8th>());
    }
    if(s=="cd6")
    {
        myDifferenters.push_back(
                                 std::make_shared<centraldifference6th>());
    }
    else if(s=="cd2")
    {
        myDifferenters.push_back(
                                 std::make_shared<centraldifference2nd>());
    }
    else
    {
        myDifferenters.push_back(
                                 std::make_shared<centraldifference6th>());
    }
    myDifferentersBND.push_back(std::make_shared<boundaryScheme>());
    
}

void nuc3d::fieldOperator3d::setTimeMethod()
{
    std::string s=MethodMap["scheme_time"];
    
    if(s=="RK3")
        myIntegrators.reset(new tvdrk3rd);
    else
        myIntegrators.reset(new tvdrk3rd);
}

void nuc3d::fieldOperator3d::reconstructionInner(
                                                 const Field& fieldIN,
                                                 int direction,
                                                 int upwind,
                                                 Field & fieldOUT)
{
    switch(direction)
    {
        case 0:
            myInteroplators[0]->interpolationInner(fieldIN,upwind,fieldOUT,bufferSize_xi);
            break;
        case 1:
            myInteroplators[1]->interpolationInner(fieldIN,upwind,fieldOUT,bufferSize_eta);
            break;
        case 2:
            myInteroplators[2]->interpolationInner(fieldIN,upwind,fieldOUT,bufferSize_zeta);
            break;
    }
}

void nuc3d::fieldOperator3d::reconstructionBoundary(
                                                    const Field &fieldIN,
                                                    const Field &boundaryL,
                                                    const Field &boundaryR,
                                                    const int direction,
                                                    const int upwind,
                                                    Field &fieldOUT,
                                                    int typeL,
                                                    int typeR)
{
    
    if((-1)==typeL)
        reconstructionBoundaryExteriorL(fieldIN,boundaryL,direction,upwind,fieldOUT);
    else
        reconstructionBoundaryInnerL(fieldIN,boundaryL,direction,upwind,fieldOUT);
    
    
    if((-1)==typeR)
        reconstructionBoundaryExteriorR(fieldIN,boundaryR,direction,upwind,fieldOUT);
    else
        reconstructionBoundaryInnerR(fieldIN,boundaryR,direction,upwind,fieldOUT);
    
}

void nuc3d::fieldOperator3d::reconstructionBoundaryInnerL(
                                                          const Field &fieldIN,
                                                          const Field &boundaryL,
                                                          const int direction,
                                                          const int upwind,
                                                          Field &fieldOUT)
{
    switch(direction)
    {
        case 0:
            myInteroplators[0]->interpolationBoundaryL(fieldIN,boundaryL,upwind,fieldOUT,bufferSize_xi);
            break;
        case 1:
            myInteroplators[1]->interpolationBoundaryL(fieldIN,boundaryL,upwind,fieldOUT,bufferSize_eta);
            break;
        case 2:
            myInteroplators[2]->interpolationBoundaryL(fieldIN,boundaryL,upwind,fieldOUT,bufferSize_zeta);
            break;
    }
    
}


void nuc3d::fieldOperator3d::reconstructionBoundaryInnerR(
                                                          const Field &fieldIN,
                                                          const Field &boundaryR,
                                                          const int direction,
                                                          const int upwind,
                                                          Field &fieldOUT)
{
    switch(direction)
    {
        case 0:
            myInteroplators[0]->interpolationBoundaryR(fieldIN,boundaryR,upwind,fieldOUT,bufferSize_xi);
            break;
        case 1:
            myInteroplators[1]->interpolationBoundaryR(fieldIN,boundaryR,upwind,fieldOUT,bufferSize_eta);
            break;
        case 2:
            myInteroplators[2]->interpolationBoundaryR(fieldIN,boundaryR,upwind,fieldOUT,bufferSize_zeta);
            break;
    }
    
}

void nuc3d::fieldOperator3d::reconstructionBoundaryExteriorL(
                                                          const Field &fieldIN,
                                                          const Field &boundaryL,
                                                          const int direction,
                                                          const int upwind,
                                                          Field &fieldOUT)
{
    switch(direction)
    {
        case 0:
            myInteroplatorsBND[0]->interpolationBoundaryL(fieldIN,boundaryL,upwind,fieldOUT,bufferSize_xi);
            break;
        case 1:
            myInteroplatorsBND[1]->interpolationBoundaryL(fieldIN,boundaryL,upwind,fieldOUT,bufferSize_eta);
            break;
        case 2:
            myInteroplatorsBND[2]->interpolationBoundaryL(fieldIN,boundaryL,upwind,fieldOUT,bufferSize_zeta);
            break;
    }
    
}


void nuc3d::fieldOperator3d::reconstructionBoundaryExteriorR(
                                                          const Field &fieldIN,
                                                          const Field &boundaryR,
                                                          const int direction,
                                                          const int upwind,
                                                          Field &fieldOUT)
{
    switch(direction)
    {
        case 0:
            myInteroplatorsBND[0]->interpolationBoundaryR(fieldIN,boundaryR,upwind,fieldOUT,bufferSize_xi);
            break;
        case 1:
            myInteroplatorsBND[1]->interpolationBoundaryR(fieldIN,boundaryR,upwind,fieldOUT,bufferSize_eta);
            break;
        case 2:
            myInteroplatorsBND[2]->interpolationBoundaryR(fieldIN,boundaryR,upwind,fieldOUT,bufferSize_zeta);
            break;
    }
    
}

void nuc3d::fieldOperator3d::differenceInner(const Field& fieldIN,
                                             int direction,
                                             Field & fieldOUT)
{
    switch(direction)
    {
        case 0:
            myDifferenters[0]->differentialInner(fieldIN,fieldOUT,bufferSize_xi);
            break;
        case 1:
            myDifferenters[1]->differentialInner(fieldIN,fieldOUT,bufferSize_eta);
            break;
        case 2:
            myDifferenters[2]->differentialInner(fieldIN,fieldOUT,bufferSize_zeta);
            break;
    }
}

void nuc3d::fieldOperator3d::differenceBoundary(
                                                const Field &fieldIN,
                                                const Field &boundaryL,
                                                const Field &boundaryR,
                                                int direction,
                                                Field &fieldOUT,
                                                int typeL,
                                                int typeR)
{
    
    if((-1)==typeL)
        differenceBoundaryExteriorL(fieldIN,boundaryL,direction,fieldOUT);
    else
        differenceBoundaryInnerL(fieldIN,boundaryL,direction,fieldOUT);
    
    if((-1)==typeR)
        differenceBoundaryExteriorR(fieldIN,boundaryR,direction,fieldOUT);
    else
        differenceBoundaryInnerR(fieldIN,boundaryR,direction,fieldOUT);
}


void nuc3d::fieldOperator3d::differenceBoundaryInnerL(const Field &fieldIN,
                                                      const Field &boundaryL,
                                                      const int direction,
                                                      Field &fieldOUT)
{
    switch(direction)
    {
        case 0:
            myDifferenters[0]->differentialBoundaryL(fieldIN,boundaryL,fieldOUT,bufferSize_xi);
            break;
        case 1:
            myDifferenters[1]->differentialBoundaryL(fieldIN,boundaryL,fieldOUT,bufferSize_eta);
            break;
        case 2:
            myDifferenters[2]->differentialBoundaryL(fieldIN,boundaryL,fieldOUT,bufferSize_zeta);
            break;
    }
}

void nuc3d::fieldOperator3d::differenceBoundaryExteriorL(const Field &fieldIN,
                                                         const Field &boundaryL,
                                                         const int direction,
                                                         Field &fieldOUT)
{
    switch(direction)
    {
        case 0:
            myDifferentersBND[0]->differentialBoundaryL(fieldIN,boundaryL,fieldOUT,bufferSize_xi);
            break;
        case 1:
            myDifferentersBND[1]->differentialBoundaryL(fieldIN,boundaryL,fieldOUT,bufferSize_eta);
            break;
        case 2:
            myDifferentersBND[2]->differentialBoundaryL(fieldIN,boundaryL,fieldOUT,bufferSize_zeta);
            break;
    }
    
    
}

void nuc3d::fieldOperator3d::differenceBoundaryInnerR(const Field &fieldIN,
                                                      const Field &boundaryR,
                                                      const int direction,
                                                      Field &fieldOUT)
{
    switch(direction)
    {
        case 0:
            myDifferenters[0]->differentialBoundaryR(fieldIN,boundaryR,fieldOUT,bufferSize_xi);
            break;
        case 1:
            myDifferenters[1]->differentialBoundaryR(fieldIN,boundaryR,fieldOUT,bufferSize_eta);
            break;
        case 2:
            myDifferenters[2]->differentialBoundaryR(fieldIN,boundaryR,fieldOUT,bufferSize_zeta);
            break;
    }
    
}

void nuc3d::fieldOperator3d::differenceBoundaryExteriorR(const Field &fieldIN,
                                                         const Field &boundaryR,
                                                         const int direction,
                                                         Field &fieldOUT)
{
    switch(direction)
    {
        case 0:
            myDifferentersBND[0]->differentialBoundaryR(fieldIN,boundaryR,fieldOUT,bufferSize_xi);
            break;
        case 1:
            myDifferentersBND[1]->differentialBoundaryR(fieldIN,boundaryR,fieldOUT,bufferSize_eta);
            break;
        case 2:
            myDifferentersBND[2]->differentialBoundaryR(fieldIN,boundaryR,fieldOUT,bufferSize_zeta);
            break;
    }
    
    
}

void nuc3d::fieldOperator3d::timeIntegral (const VectorField &rhs, // rhs
                                           const VectorField &u_n, // u_n
                                           VectorField &u_i, //u_i
                                           double dt,
                                           int step)
{
    myIntegrators->integrationAll(rhs,u_n,u_i,dt,step);
}

#endif
