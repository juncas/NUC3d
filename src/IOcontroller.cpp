#ifndef IOController_cpp
#define IOController_cpp
#include "IOController.h"

/**************************************************************************************
Definition of class IOController: 
			This class is used to store numerical method input parameters.
**************************************************************************************/
/**************************************************************************************
						Definition of constructors and destructors
**************************************************************************************/
nuc3d::IOController::IOController():
    myIOController
    (
     {
         {"startStep",0},
         {"endStep",100000},
         {"saveStep",1000},
         {"currentStep",0}
     }
    ),
    myTimeController
    (
     {
         {"cfl",0.5},
         {"tl",0.0},
         {"dt",0.01},
         {"currentTime",0.0}
     }
    )

{
	std::ifstream file("IOController.in");

	auto count=myTimeController.size()+myIOController.size();	if(file)
	{
		while(readIOFile(file,myIOController,myTimeController))
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
};

nuc3d::IOController::~IOController()
{

};
/**************************************************************************************
						Definition of member functions
**************************************************************************************/
std::istream& nuc3d::IOController::readIOFile(std::istream& ios, std::map<std::string,int> &Methods, std::map<std::string,double> &Time)
{
	std::string word0;
	std::string word1;
	ios>>word0>>word1;

	if(Methods.find(word0)!=Methods.end())
	{
		std::istringstream value0(word1);
		int value;
		value0>>value;
		Methods[word0]=value;
	}
	else if(Time.find(word0)!=Time.end())
	{
		std::istringstream value0(word1);
		double value;
		value0>>value;
		Time[word0]=value;
	}
	else if(word0.size())
	{
		std::cout<<"Parameter \'"<<word0<<"in file \'IOController.io\' is not supported!"<<std::endl;
		std::cout<<"Make sure that every parameter exist in the following name list:"<<std::endl;
		for(auto beg=Methods.begin();beg!=Methods.end();beg++)
			std::cout<<beg->first<<" "<<beg->second<<std::endl;
		for(auto beg=Time.begin();beg!=Time.end();beg++)
			std::cout<<beg->first<<" "<<beg->second<<std::endl;
	}
	return ios;
}


double nuc3d::IOController::getValue(std::string s)
{	
	if(myTimeController.find(s)!=myTimeController.end())
	{
		return myTimeController[s];	
	}
	else
	{		
		std::cout<<"ERROR:Parameter \'"<<s<<" you are requring does not exist!"<<std::endl;
		exit(0);
	}
};



int nuc3d::IOController::getStep(std::string s)
{
    if(myIOController.find(s)!=myIOController.end())
    {
        return myIOController[s];
    }
    else
    {
        std::cout<<"ERROR:Parameter \'"<<s<<" you are requring does not exist!"<<std::endl;
        exit(0);
    }
};

bool nuc3d::IOController::ifsolve()
{
    double currentTime=myTimeController["currentTime"];
    double tl=myTimeController["tl"];
    
    int currentStep=myIOController["currentStep"];
    int endStep=myIOController["endStep"];
    
    if((currentTime<=tl)&&(currentStep<endStep))
        return true;
    else
        return false;
}

bool nuc3d::IOController::ifsave()
{
    int saveStep=myIOController["saveStep"];
    int currentStep=myIOController["currentStep"];
    
    if(currentStep%saveStep)
        return true;
    else
        return false;
}
/**************************************************************************************
								End of definition
**************************************************************************************/
#endif
