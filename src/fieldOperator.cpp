#ifndef fieldOperator3d_cpp
#define fieldOperator3d_cpp
#include "fieldOperator.h"
#include "weno5js.h"
#include "TVD-RK3.h"
#include "centraldifference6th.h"

nuc3d::fieldOperator3d::fieldOperator3d(const VectorField &U):
MethodMap
(
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
    std::ifstream file("IOController.in");
    auto count=MethodMap.size();
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
    
    setTimeMethod(U);
    
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
    
};

nuc3d::fieldOperator3d::~fieldOperator3d()
{
    
};
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
        std::cout<<"Make sure that every parameter exist in the following name list:"<<std::endl;
        for(auto beg=Methods.begin();beg!=Methods.end();beg++)
            std::cout<<beg->first<<" "<<beg->second<<std::endl;
    }
    return ios;
}


void nuc3d::fieldOperator3d::setMethodIvsX()
{
    std::string s=MethodMap["scheme_x_ivs"];
    
    if(s=="weno5js")
    {
        bufferSize=3;
        myInteroplators.push_back(std::make_shared<interoplation>(new weno5js));
    }
    else
        myInteroplators.push_back(std::make_shared<interoplation>(new weno5js));
}

void nuc3d::fieldOperator3d::setMethodIvsY()
{
    std::string s=MethodMap["scheme_y_ivs"];
    
    if(s=="weno5js")
    {
        bufferSize=3;
        myInteroplators.push_back(std::make_shared<interoplation>(new weno5js));
    }
    else
        myInteroplators.push_back(std::make_shared<interoplation>(new weno5js));}

void nuc3d::fieldOperator3d::setMethodIvsZ()
{
    std::string s=MethodMap["scheme_z_ivs"];
    
    if(s=="weno5js")
        myInteroplators.push_back(std::make_shared<interoplation>(new weno5js));
    else
        myInteroplators.push_back(std::make_shared<interoplation>(new weno5js));
}

void nuc3d::fieldOperator3d::setDiffMethodX()
{
    std::string s=MethodMap["scheme_x_vis"];
    
    if(s=="cd6")
    {
        bufferSize=3;
        myDifferenters.push_back(
                                 std::make_shared<differential>(new centraldifference6th));
    }
    else
        myDifferenters.push_back(
                                 std::make_shared<differential>(new centraldifference6th));}

void nuc3d::fieldOperator3d::setDiffMethodY()
{
    std::string s=MethodMap["scheme_y_vis"];
    
    if(s=="cd6")
    {
        bufferSize=3;
        myDifferenters.push_back(
                                 std::make_shared<differential>(new centraldifference6th));
    }
    else
        myDifferenters.push_back(
                                 std::make_shared<differential>(new centraldifference6th));
}

void nuc3d::fieldOperator3d::setDiffMethodZ()
{
    std::string s=MethodMap["scheme_z_vis"];
    
    if(s=="cd6")
    {
        bufferSize=3;
        
        myDifferenters.push_back(
                                 std::make_shared<differential>(new centraldifference6th));
    }
    else
        myDifferenters.push_back(
                                 std::make_shared<differential>(new centraldifference6th));
}

void nuc3d::fieldOperator3d::setTimeMethod(const VectorField &u)
{
    std::string s=MethodMap["scheme_time"];
    
    if(s=="RK3")
        std::shared_ptr<integration> p(new tvdrk3rd);
    else
        myIntegrators=new tvdrk3(u);
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
            myInteroplators[0]->interpolationInner(fieldIN,1,0,0,upwind,fieldOUT);
            break;
        case 1:
            myInteroplators[1]->interpolationInner(fieldIN,0,1,0,upwind,fieldOUT);
            break;
        case 2:
            myInteroplators[2]->interpolationInner(fieldIN,0,0,1,upwind,fieldOUT);
            break;
    }
}

void nuc3d::fieldOperator3d::reconstructionBoundary(
                                                    const Field &fieldIN,
                                                    const Field &boundaryL,
                                                    const Field &boundaryR,
                                                    const int direction,
                                                    const int upwind,
                                                    Field &fieldOUT)
{
    switch(direction)
    {
        case 0:
            myInteroplators[0]->interpolationInner(fieldIN,boundaryL,boundaryR,1,0,0,upwind,fieldOUT);
            break;
        case 1:
            myInteroplators[1]->interpolationInner(fieldIN,boundaryL,boundaryR,0,1,0,upwind,fieldOUT);
            break;
        case 2:
            myInteroplators[2]->interpolationInner(fieldIN,boundaryL,boundaryR,0,0,1,upwind,fieldOUT);
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
            myDifferenters[0]->interpolationInner(fieldIN,1,0,0,fieldOUT);
            break;
        case 1:
            myDifferenters[1]->interpolationInner(fieldIN,0,1,0,fieldOUT);
            break;
        case 2:
            myDifferenters[2]->interpolationInner(fieldIN,0,0,1,fieldOUT);
            break;
    }
}

void nuc3d::fieldOperator3d::differenceBoundary(
                                                const Field &fieldIN,
                                                const Field &boundaryL,
                                                const Field &boundaryR,
                                                const int direction,
                                                Field &fieldOUT)
{
    switch(direction)
    {
        case 0:
            myDifferenters[0]->interpolationInner(fieldIN,boundaryL,boundaryR,1,0,0,fieldOUT);
            break;
        case 1:
            myDifferenters[1]->interpolationInner(fieldIN,boundaryL,boundaryR,0,1,0,fieldOUT);
            break;
        case 2:
            myDifferenters[2]->interpolationInner(fieldIN,boundaryL,boundaryR,0,0,1,fieldOUT);
            break;
    }
    
}


void nuc3d::fieldOperator3d::timeIntegral (
                                           const Field& dfdx, //dfdx
                                           const Field& dfdy, //dfdy
                                           const Field& dfdz, //dfdz
                                           const Field& u0, // u0
                                           const Field& un, // un
                                           const Field& rhs, // rhs
                                           const double dt,
                                           Field &fieldOUT,
                                           int step)
{
    myIntegrators->getRHS(dfdx,dfdy,dfdz,dt, rhs);
    myIntegrators->integrationAll(rhs,u0,un,step,fieldOUT);
    
}
#endif
