#ifndef cell_h
#define cell_h

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <vector>

namespace nuc3d
{
	class cell
	{
		std::vector<double> data;
	public:
		cell(int);
		~cell();
	};
}

#endif