#ifndef IOController_h
#define IOController_h

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
namespace nuc3d
{
	class IOController
	{
	public:
		double cfl;
		double tl;
		
        std::map<std::string,int> myIOController;
        std::map<std::string,double> myTimeController;
	public:
        
		IOController();
		~IOController();
        bool ifsolve();
        bool ifsave();
        void increaseStep();
        
	public:
		int getMethod(std::string);
		double getValue(std::string);
		std::istream& readIOFile(std::istream&,std::map<std::string,int> &,std::map<std::string,double> &);
        
        
	};
}
#endif