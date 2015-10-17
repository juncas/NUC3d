#ifndef meshGrid_cpp
#define meshGrid_cpp
#include "meshgrid.h" 
#include <sstream>


/**************************************************************************************
			Definition of class communicator: base class for communicators
**************************************************************************************/
/**************************************************************************************
						Definition of constructors and destructors
**************************************************************************************/
nuc3d::meshGrid::meshGrid(int& argc, char **& argv) :
	myComm(argc, argv), myCtrler(), myPhys(), nblocks(0), myOperator()
{
	input();
	initial();
}
nuc3d::meshGrid::~meshGrid()
{

}
/**************************************************************************************
						Definition of member functions
**************************************************************************************/
void nuc3d::meshGrid::input()
{
	std::stringstream ss;
	int ProcId = myComm.getMyId();
	int nx0, ny0, nz0;
	std::string forename = ("./meshes/mesh_");
	std::string midname;
	std::string tailname = (".dat");

	ss << ProcId;
	ss >> midname;

	std::string filename = forename + midname + tailname;


	std::ifstream myFile(filename);
	if (myFile)
	{
		while (myFile >> nx0 >> ny0 >> nz0)
		{
			myBlocks.push_back(MeshBlock(nx0, ny0, nz0, myOperator.getBufferSize(), myPhys.getEqNum()));

			if (!myBlocks.empty())
			{
				auto block = myBlocks.back();
				block.setMyId(nblocks);
				block.input(myFile);
				nblocks++;
			}
			else
			{
				std::cout << nblocks << " blocks have been read! Failed" << std::endl;
				myComm.AbortMPI();
			}
		};

		std::cout << nblocks << " blocks have been read!" << std::endl;
	}
	else
	{
		std::cout << "file " << filename << " does not exist" << std::endl;
		myComm.AbortMPI();
	}
	myFile.close();
}

void nuc3d::meshGrid::initial()
{
	for (auto iterBlocks = myBlocks.begin(); iterBlocks != myBlocks.end(); ++iterBlocks)
	{
		iterBlocks->initial();
		myPhys.initial(iterBlocks->myEuler);
	}
}

void nuc3d::meshGrid::solve()
{
    while (myCtrler.solve()) {
        int rkstep=0;
        
        while(myOperator.myIntegrators->nstep>rkstep)
        {
            solveRiemann();
            solveBounaryConditions();
            solveInvicidFlux();
            solveViscousFLux();
            solveGetRHS();
            solveIntegral();
        }
        
    }
    
}

void nuc3d::meshGrid::solveRiemann()
{
	for (auto iterBlocks = myBlocks.begin(); iterBlocks != myBlocks.end(); ++iterBlocks)
	{
		myPhys.solve(iterBlocks->myEuler);
	}

}

void nuc3d::meshGrid::solveInvicidFlux()
{

	for (int i = 0; i < myPhys.getEqNum(); i++)
	{
		for (auto iterBlocks = myBlocks.begin(); iterBlocks != myBlocks.end(); ++iterBlocks)
		{
			auto pEuler = iterBlocks->myEuler;
			iterBlocks->myBuffers[i].setBufferSend(pEuler.FluxL_xi[i]);
			myComm.bufferSendRecv(iterBlocks->myBuffers[i], iterBlocks - myBlocks.begin());
			myOperator.reconstructionInner(pEuler.FluxL_xi[i], 0, 1, pEuler.reconstFluxL_xi[i]);

		}
		for (auto iterBlocks = myBlocks.begin(); iterBlocks != myBlocks.end(); ++iterBlocks)
		{
			auto pEuler = iterBlocks->myEuler;
			myComm.waitAllSendRecv(iterBlocks->myBuffers[i]);
			myOperator.reconstructionBoundary(pEuler.FluxL_xi[i],
				(iterBlocks->myBuffers[i].BufferRecv[0]),
				(iterBlocks->myBuffers[i].BufferRecv[1]),
				0,
				1,
				pEuler.reconstFluxL_xi[i]);
		}

		for (auto iterBlocks = myBlocks.begin(); iterBlocks != myBlocks.end(); ++iterBlocks)
		{
			auto pEuler = iterBlocks->myEuler;
			iterBlocks->myBuffers[i].setBufferSend(pEuler.FluxL_eta[i]);
			myComm.bufferSendRecv(iterBlocks->myBuffers[i], iterBlocks - myBlocks.begin());
			myOperator.reconstructionInner(pEuler.FluxL_eta[i], 1, 1, pEuler.reconstFluxL_eta[i]);

		}
		for (auto iterBlocks = myBlocks.begin(); iterBlocks != myBlocks.end(); ++iterBlocks)
		{
			auto pEuler = iterBlocks->myEuler;
			myComm.waitAllSendRecv(iterBlocks->myBuffers[i]);
			myOperator.reconstructionBoundary(pEuler.FluxL_eta[i],
				(iterBlocks->myBuffers[i].BufferRecv[2]),
				(iterBlocks->myBuffers[i].BufferRecv[3]),
				1,
				1,
				pEuler.reconstFluxL_eta[i]);
		}

		for (auto iterBlocks = myBlocks.begin(); iterBlocks != myBlocks.end(); ++iterBlocks)
		{
			auto pEuler = iterBlocks->myEuler;
			iterBlocks->myBuffers[i].setBufferSend(pEuler.FluxL_zeta[i]);
			myComm.bufferSendRecv(iterBlocks->myBuffers[i], iterBlocks - myBlocks.begin());
			myOperator.reconstructionInner(pEuler.FluxL_eta[i], 2, 1, pEuler.reconstFluxL_eta[i]);

		}
		for (auto iterBlocks = myBlocks.begin(); iterBlocks != myBlocks.end(); ++iterBlocks)
		{
			auto pEuler = iterBlocks->myEuler;
			myComm.waitAllSendRecv(iterBlocks->myBuffers[i]);
			myOperator.reconstructionBoundary(pEuler.FluxL_eta[i],
				(iterBlocks->myBuffers[i].BufferRecv[4]),
				(iterBlocks->myBuffers[i].BufferRecv[5]),
				2,
				1,
				pEuler.reconstFluxL_eta[i]);
		}
	}

	for (int i = 0; i < myPhys.getEqNum(); i++)
	{
		for (auto iterBlocks = myBlocks.begin(); iterBlocks != myBlocks.end(); ++iterBlocks)
		{
			auto pEuler = iterBlocks->myEuler;
			iterBlocks->myBuffers[i].setBufferSend(pEuler.FluxR_xi[i]);
			myComm.bufferSendRecv(iterBlocks->myBuffers[i], iterBlocks - myBlocks.begin());
			myOperator.reconstructionInner(pEuler.FluxR_xi[i], 0, -1, pEuler.reconstFluxR_xi[i]);

		}
		for (auto iterBlocks = myBlocks.begin(); iterBlocks != myBlocks.end(); ++iterBlocks)
		{
			auto pEuler = iterBlocks->myEuler;
			myComm.waitAllSendRecv(iterBlocks->myBuffers[i]);
			myOperator.reconstructionBoundary(pEuler.FluxR_xi[i],
				(iterBlocks->myBuffers[i].BufferRecv[0]),
				(iterBlocks->myBuffers[i].BufferRecv[1]),
				0,
				-1,
				pEuler.reconstFluxR_xi[i]);
		}

		for (auto iterBlocks = myBlocks.begin(); iterBlocks != myBlocks.end(); ++iterBlocks)
		{
			auto pEuler = iterBlocks->myEuler;
			iterBlocks->myBuffers[i].setBufferSend(pEuler.FluxR_eta[i]);
			myComm.bufferSendRecv(iterBlocks->myBuffers[i], iterBlocks - myBlocks.begin());
			myOperator.reconstructionInner(pEuler.FluxR_eta[i], 1, -1, pEuler.reconstFluxR_eta[i]);

		}
		for (auto iterBlocks = myBlocks.begin(); iterBlocks != myBlocks.end(); ++iterBlocks)
		{
			auto pEuler = iterBlocks->myEuler;
			myComm.waitAllSendRecv(iterBlocks->myBuffers[i]);
			myOperator.reconstructionBoundary(pEuler.FluxR_eta[i],
				(iterBlocks->myBuffers[i].BufferRecv[2]),
				(iterBlocks->myBuffers[i].BufferRecv[3]),
				1,
				-1,
				pEuler.reconstFluxR_eta[i]);
		}

		for (auto iterBlocks = myBlocks.begin(); iterBlocks != myBlocks.end(); ++iterBlocks)
		{
			auto pEuler = iterBlocks->myEuler;
			iterBlocks->myBuffers[i].setBufferSend(pEuler.FluxR_zeta[i]);
			myComm.bufferSendRecv(iterBlocks->myBuffers[i], iterBlocks - myBlocks.begin());
			myOperator.reconstructionInner(pEuler.FluxR_zeta[i], 2, -1, pEuler.reconstFluxR_zeta[i]);

		}
		for (auto iterBlocks = myBlocks.begin(); iterBlocks != myBlocks.end(); ++iterBlocks)
		{
			auto pEuler = iterBlocks->myEuler;
			myComm.waitAllSendRecv(iterBlocks->myBuffers[i]);
			myOperator.reconstructionBoundary(pEuler.FluxR_zeta[i],
				(iterBlocks->myBuffers[i].BufferRecv[4]),
				(iterBlocks->myBuffers[i].BufferRecv[5]),
				2,
				-1,
				pEuler.reconstFluxR_zeta[i]);
		}
	}

}

void nuc3d::meshGrid::solveViscousFLux()
{
    
}

void nuc3d::meshGrid::solveGetRHS()
{
    
}
/**************************************************************************************
								End of definition
**************************************************************************************/

#endif