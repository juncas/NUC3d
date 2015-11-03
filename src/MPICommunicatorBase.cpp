#ifndef MPICommunicator_cpp
#define MPICommunicator_cpp
#include "MPICommunicatorBase.h"
/**************************************************************************************
			Definition of class communicator: base class for communicators
**************************************************************************************/
/**************************************************************************************
						Definition of constructors and destructors
**************************************************************************************/
nuc3d::MPICommunicator::MPICommunicator()
{
    MPI_Comm_rank(MPI_COMM_WORLD, &myProc);
  	MPI_Comm_size(MPI_COMM_WORLD, &nProc);
}

nuc3d::MPICommunicator::~MPICommunicator()
{
  MPI_Finalize ();
}
/**************************************************************************************
						Definition of member functions
**************************************************************************************/
void nuc3d::MPICommunicator::FinializeMPI()
{
  MPI_Finalize ();
}
void nuc3d::MPICommunicator::AbortMPI()
{
  	MPI_Abort ( MPI_COMM_WORLD,99 );
}
/**************************************************************************************
								End of definition
**************************************************************************************/
#endif