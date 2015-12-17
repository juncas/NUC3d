#ifndef field_h
#define field_h

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <vector>
namespace nuc3d
{
	template<typename T> class field3d
	{
	private:
		T *field_data;

		int nx;
		int ny;
		int nz;

	public:
		//Constructor 

			//copy constructor
		field3d(const field3d<T> &);

		field3d(int, int, int);

		field3d(int, int, int, T);

		//Destructor 
		~field3d();

	public://Member functions
		int getSizeX() const
		{
			return nx;
		};
		int getSizeY() const
		{
			return ny;
		};
		int getSizeZ() const
		{
			return nz;
		};
		int getSize() const
		{
			return nx*ny*nz;
		}
		T getValue(int, int, int) const;
		void setValue(int, int, int, T);
		int getIndex(int, int, int) const;
		void resetData(int, int, int);
		T* getDataPtr() const{ return field_data; } ;
		static const int fieldDim = 3;
	public://Operators
		field3d<T>& operator=(const field3d<T> &);

		field3d<T>& operator+=(const field3d<T> &);
		field3d<T>& operator-=(const field3d<T> &);

		field3d<T> operator+(const field3d<T> &);
		field3d<T> operator-(const field3d<T> &);

		field3d<T>& operator*=(const field3d<T> &);
		field3d<T> operator* (const field3d<T> &);

		field3d<T>& operator*=(const T &);
		field3d<T> operator* (const T &);


		field3d<T>& operator/=(const field3d<T> &);
		field3d<T>  operator/ (const field3d<T> &);

		field3d<T>& operator/=(const T &);
		field3d<T>  operator/ (const T &);

	};

	typedef field3d<double> Field;
	typedef std::vector<Field> VectorField;

}

#include "field.cpp"

#endif