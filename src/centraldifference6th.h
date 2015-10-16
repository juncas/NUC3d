#ifndef centraldifference6th_h
#define centraldifference6th_h
#include "fieldOperator.h"


namespace nuc3d
{
	class centraldifference6th : public differential
	{
		const double coeff_cd6_alpha[3] = { 3.0 / 4.0, -3.0 / 20.0,1.0 / 60.0 };

    public:
    	centraldifference6th();
    	~centraldifference6th();
    	 // 6th order central difference functions

        void differentialInner(
                const Field &,
                const int,
                const int,
                const int,
                Field &);
        void differentialBoundary(
                const Field &,
                const Field &,
                const Field &,
                const int,
                const int,
                const int,
                Field &);
    private:
        double cd6differential(const double *);

	};
}



#endif