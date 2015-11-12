//
//  boundaryConditions.hpp
//  NUC3d
//
//  Created by Jun Peng on 15/11/5.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#ifndef boundaryConditions_hpp
#define boundaryConditions_hpp
#include <memory>
#include <map>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include "field.h"
#include "mpi.h"
#include "bufferData.hpp"

namespace nuc3d
{
	class EulerData3D;
	class MPIComunicator3d_nonblocking;
	class PDEData3d;
	class physicsModel;
	class EulerFlux;

	class faceBC
	{
		public:
			int Type; // -1:= bc  n:= neighbour block id
			int id; // bc id or neighbour face id

			faceBC(int myType,int myID):Type(myType),id(myID){};
			~faceBC(){};
	};

	class boundaryCondition
	{

		typedef void (boundaryCondition::*psetBuffer)(VectorBuffer &,VectorField &,int);

		psetBuffer mySetter[5]={
			&boundaryCondition::setBuffer_inlet,
			&boundaryCondition::setBuffer_outlet,
			&boundaryCondition::setBuffer_wall,
			&boundaryCondition::setBuffer_symm,
			&boundaryCondition::setBuffer_period
		};

		std::vector<faceBC> BCTopo;

		std::vector<std::vector<double>> BCvalue;//(rho,rhou,rhow,rhoe)

		public:
		boundaryCondition();
		~boundaryCondition();

		public:
		void setBC(PDEData3d &,
				physicsModel &,
				EulerData3D &);

		void updateBuffer_xi(VectorBuffer &,VectorField &);
		void updateBuffer_eta(VectorBuffer &,VectorField &);
		void updateBuffer_zeta(VectorBuffer &,VectorField &);

		void initialBC(VectorBuffer &,
				MPIComunicator3d_nonblocking &);
		private:
		std::ifstream& readBCTopo(std::ifstream&,int);

		void setBC_Outlet();
		void setBC_Wall();
		void setBC_Periodic();
		void setBC_Symmetric();

		//intlet
		void setBC_Inlet(PDEData3d &,
				EulerData3D &,
				physicsModel &myMod,
				int);

		void setBuffer_inlet(VectorBuffer &,
				VectorField &,
				int);

		void setBuffer_inlet_xi(VectorBuffer &,
				VectorField &,
				int );

		void setBuffer_inlet_eta(VectorBuffer &,
				VectorField &,
				int );

		void setBuffer_inlet_zeta(VectorBuffer &,
				VectorField &,
				int );

		//outlet
		void setBC_Outlet(PDEData3d &,
				EulerData3D &,
				physicsModel &myMod,
				int);

		void setBuffer_outlet(VectorBuffer &,
				VectorField &,
				int);

		void setBuffer_outlet_xi(VectorBuffer &,
				VectorField &,
				int);

		void setBuffer_outlet_eta(VectorBuffer &,
				VectorField &,
				int);

		void setBuffer_outlet_zeta(VectorBuffer &,
				VectorField &,
				int);

		//wall
		void setBC_wall(PDEData3d &,
				EulerData3D &,
				physicsModel &myMod,
				int);

		void setBuffer_wall(VectorBuffer &,
				VectorField &,
				int);

		void setBuffer_wall_xi(VectorBuffer &,
				VectorField &,
				int);

		void setBuffer_wall_eta(VectorBuffer &,
				VectorField &,
				int);

		void setBuffer_wall_zeta(VectorBuffer &,
				VectorField &,
				int);


		//symmetric
		void setBC_symm(PDEData3d &,
				EulerData3D &,
				physicsModel &myMod,
				int);

		void setBuffer_symm(VectorBuffer &,
				VectorField &,
				int);

		void setBuffer_symm_xi(VectorBuffer &,
				VectorField &,
				int);

		void setBuffer_symm_eta(VectorBuffer &,
				VectorField &,
				int);

		void setBuffer_symm_zeta(VectorBuffer &,
				VectorField &,
				int);
		//periodic
		void setBC_period(PDEData3d &,
				EulerData3D &,
				physicsModel &myMod,
				int);

		void setBuffer_period(VectorBuffer &,
				VectorField &,
				int);

		void setBuffer_period_xi(VectorBuffer &,
				VectorField &,
				int);

		void setBuffer_period_eta(VectorBuffer &,
				VectorField &,
				int);

		void setBuffer_period_zeta(VectorBuffer &,
				VectorField &,
				int);

	};
}

#endif /* boundaryConditions_hpp */
