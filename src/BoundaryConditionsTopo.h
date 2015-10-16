#ifndef BoundaryConditionsTopo_h
#define BoundaryConditionsTopo_h

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>

namespace nuc3d
{
	struct BoundaryConditionsTopo
	{
		const std::string BCnames[] = { "Inlet",
			"Outlet",
			"Interior",
			"ISO_Thermal_wall",
			"Symmetric",
			"User_defined_BC" };

		int BoundaryConditions[6];
		int ConnectiveInfo[6];

	};
}

#endif