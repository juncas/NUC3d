#ifndef meshGrid_h
#define meshGrid_h

#include <cstdlib>
#include <iostream>
#include <vector>

#include "meshblock.h"
#include "MPICommunicator.h"
#include "IOcontroller.h"
#include "fieldOperator.h"
#include "physicsModel.h"


namespace nuc3d 
{
	class meshGrid
	{
		double dataLocal;
		int nblocks;

		MPIComunicator3d_nonblocking myComm; //finished
		IOController myCtrler;//finished
		physicsModel myPhys;//finished
		fieldOperator3d myOperator;		

		std::vector<MeshBlock> myBlocks;//In this version only one grid in each domain
								   //Field data are stored in these blocks

	public:

		meshGrid(int&,char **&);
		~meshGrid();

	public:

		void solve();

	private:
		void solveFlux();
		void solveViscous();
		void solveGetRHS();
		void solveIntegral();

		void MaxReduceMPI();
		void MaxReduceLocal();

		void MinReduceMPI();
		void MinReduceLocal();

		void input();
		void initial();
		void output();
		
		int getBlockNum() const {return myBlocks.size();};

	};
}

#endif