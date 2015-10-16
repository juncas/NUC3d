#ifndef MPICommunicator_h
#define MPICommunicator_h

# include <cstdlib>
# include <iostream>
# include "mpi.h"
namespace nuc3d 
{
	class MPICommunicator// resposible for transportation only
	{
		int myProc;
		int nProc;		
	public:
		MPICommunicator(int&,char **&);
		~MPICommunicator();
	public:
		int getMyId() const {return myProc;};
		int getSize() const {return nProc;};
		void FinializeMPI();
		void AbortMPI();
	};
}
#endif