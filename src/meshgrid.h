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

		MPIComunicator3d_nonblocking myComm;
		IOController myCtrler;
		physicsModel myPhys;
		std::vector<MeshBlock> myBlocks;

	public:

		meshGrid(int&,char **&);
		~meshGrid();

	public:

		void solve();

	private:
        void solveRiemann();
		void solveInvicidFlux();
		void solveViscousFLux();
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