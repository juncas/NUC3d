#ifndef field_CPP
#define field_CPP
#include "field.h"
/**************************************************************************************
						Definition of constructors and destructors
**************************************************************************************/
template<typename T> 
nuc3d::field3d<T>::field3d (const field3d<T> &field0)
{
	auto i0=field0.getSizeX();
	auto j0=field0.getSizeY();
	auto k0=field0.getSizeZ();

	nx=i0;
	ny=j0;
	nz=k0;

	field_data=new T[nx*ny*nz];

	if (!(field_data==nullptr))
	{
		for(int k=0;k!=nz;++k)
			for(int j=0;j!=ny;++j)
				for(int i=0;i!=nx;++i)
					{
						int idx_loc=getIndex(i,j,k);
						field_data[idx_loc]=field0.field_data[idx_loc];
					}
	}
	else
	{
		std::cout<<"field memory allocation failed !\n"<<std::endl;
	}

};

template<typename T> 
nuc3d::field3d<T>::field3d (int i0,int j0,int k0):
	nx(i0),ny(j0),nz(k0) 
{
	field_data=new T[i0*j0*k0];

	if (!(field_data==nullptr))
	{
		for(int k=0;k!=nz;++k)
			for(int j=0;j!=ny;++j)
				for(int i=0;i!=nx;++i)
					{
						int idx_loc=getIndex(i,j,k);

						//this line is unsafe!!
						field_data[idx_loc]=1.0;
					}
	}
	else
	{
		std::cout<<"field memory allocation failed !\n"<<std::endl;
	}

};

template<typename T> 
nuc3d::field3d<T>::field3d (int i0,int j0,int k0,T val):
nx(i0),ny(j0),nz(k0) 
{		
	field_data=new T[i0*j0*k0];
	if (!(field_data==nullptr))
	{
		for(int k=0;k!=nz;++k)
			for(int j=0;j!=ny;++j)
				for(int i=0;i!=nx;++i)
					{
						int idx_loc=getIndex(i,j,k);
						field_data[idx_loc]=val;
					}
	}
	else
	{
		std::cout<<"field memory allocation failed !\n"<<std::endl;
	}

};

template<typename T> 
nuc3d::field3d<T>::~field3d()
{
	if(!(field_data==nullptr))
		delete [] field_data;
};


/**************************************************************************************
						Definition of member functions
**************************************************************************************/
template<typename T> 
int nuc3d::field3d<T>::getIndex(int i,int j,int k) const
{
	return k*nx*ny+j*nx+i;
};

template<typename T> 
T nuc3d::field3d<T>::getValue(int i,int j,int k) const
{
	return field_data[getIndex(i,j,k)];
};

template<typename T> 
void nuc3d::field3d<T>::setValue(int i,int j,int k,T val)
{
	field_data[getIndex(i,j,k)]=val;
};

template<typename T> 
void nuc3d::field3d<T>::resetData(int i,int j,int k)
{
	if(field_data==nullptr)
	{
		nx=i;
		ny=j;
		nz=k;
		field_data=new T[nx*ny*nz];
	}
	else
	{
		std::cout<<"WARNNING: Trying to reset field data!"<<std::endl;
		delete [] field_data;
		nx=i;
		ny=j;
		nz=k;
		field_data=new T[nx*ny*nz];
	}
};

/**************************************************************************************
						Definition of = operator: field A = field B
**************************************************************************************/
template<typename T> 
nuc3d::field3d<T>& nuc3d::field3d<T>::operator= (const field3d<T> &rhs)
{
	if (nx==rhs.getSizeX()&&
		ny==rhs.getSizeY()&&
		nz==rhs.getSizeZ())
	{
		if (!(field_data==nullptr)&&!(rhs.field_data==nullptr))
		{
			for(int k=0;k!=nz;++k)
				for(int j=0;j!=ny;++j)
					for(int i=0;i!=nx;++i)
					{
						int idx_loc=getIndex(i,j,k);
						field_data[idx_loc]=rhs.field_data[idx_loc];
					}
		}
		else
		{
			std::cout<<"ERROR: \n"<<std::endl;
			std::cout<<"Null field for operation '=' !\n"<<std::endl;		
		}
	}
	else
	{
		std::cout<<"ERROR: \n"<<std::endl;
		std::cout<<"Field size does not match for operation '=' !\n"<<std::endl;		
	}

	return *this;
};

/**************************************************************************************
			Definition of +,- operatiors : field A + field B,field A - field B
**************************************************************************************/
/**************************************************************************************
					1. Definition of + operatiors : field A + field B
**************************************************************************************/
template<typename T> 
nuc3d::field3d<T>& nuc3d::field3d<T>::operator+= (const field3d<T> &rhs)
{
	if (nx==rhs.getSizeX()&&
		ny==rhs.getSizeY()&&
		nz==rhs.getSizeZ())
	{
		if (!(field_data==nullptr)&&!(rhs.field_data==nullptr))
		{
			for(int k=0;k!=nz;++k)
				for(int j=0;j!=ny;++j)
					for(int i=0;i!=nx;++i)
					{
						int idx_loc=getIndex(i,j,k);
						field_data[idx_loc]+=rhs.field_data[idx_loc];
					}
		}
		else
		{
			std::cout<<"ERROR: \n"<<std::endl;
			std::cout<<"Null field for operation '+=' !\n"<<std::endl;		
		}
	}
	else
	{
		std::cout<<"ERROR: \n"<<std::endl;
		std::cout<<"Field size does not match for operation '+=' !\n"<<std::endl;		
	}

	return *this;
};

template<typename T> 
nuc3d::field3d<T> nuc3d::field3d<T>::operator+ (const field3d<T> &rhs)
{
	field3d tmp(*this);
	tmp+=rhs;

	return tmp;
};


/**************************************************************************************
					2. Definition of - operatiors : field A - field B
**************************************************************************************/
template<typename T> 
nuc3d::field3d<T>& nuc3d::field3d<T>::operator-= (const field3d<T> &rhs)
{
	if (nx==rhs.getSizeX()&&
		ny==rhs.getSizeY()&&
		nz==rhs.getSizeZ())
	{
		if (!(field_data==nullptr)&&!(rhs.field_data==nullptr))
		{
			for(int k=0;k!=nz;++k)
				for(int j=0;j!=ny;++j)
					for(int i=0;i!=nx;++i)
					{
						int idx_loc=getIndex(i,j,k);
						field_data[idx_loc]-=rhs.field_data[idx_loc];
					}
		}
		else
		{
			std::cout<<"ERROR: \n"<<std::endl;
			std::cout<<"Null field for operation '-=' !\n"<<std::endl;		
		}
	}
	else
	{
		std::cout<<"ERROR: \n"<<std::endl;
		std::cout<<"Field size does not match for operation '-=' !\n"<<std::endl;		
	}

	return *this;
};

template<typename T> 
nuc3d::field3d<T> nuc3d::field3d<T>::operator- (const field3d<T> &rhs)
{
	field3d tmp(*this);
	tmp-=rhs;

	return tmp;
};

/**************************************************************************************
							Definition of *,/ operators
**************************************************************************************/

/**************************************************************************************
				1. Definition of * operators (field,constance),(field,field)
**************************************************************************************/
/**************************************************************************************
				1.1 Definition of * operators (field,field)
**************************************************************************************/
template<typename T> 
nuc3d::field3d<T>& nuc3d::field3d<T>::operator*= (const field3d<T> &rhs)
{
	if (nx==rhs.getSizeX()&&
		ny==rhs.getSizeY()&&
		nz==rhs.getSizeZ())
	{
		if (!(field_data==nullptr)&&!(rhs.field_data==nullptr))
		{
			for(int k=0;k!=nz;++k)
				for(int j=0;j!=ny;++j)
					for(int i=0;i!=nx;++i)
					{
						int idx_loc=getIndex(i,j,k);
						field_data[idx_loc]*=rhs.field_data[idx_loc];
					}
		}
		else
		{
			std::cout<<"ERROR: \n"<<std::endl;
			std::cout<<"Null field for operation '*=' !\n"<<std::endl;		
		}
	}
	else
	{
		std::cout<<"ERROR: \n"<<std::endl;
		std::cout<<"Field size does not match for operation '*=' !\n"<<std::endl;		
	}

	return *this;
};

template<typename T> 
nuc3d::field3d<T> nuc3d::field3d<T>::operator* (const field3d<T> &rhs)
{
	field3d tmp(*this);
	tmp*=rhs;

	return tmp;
};

/**************************************************************************************
				1.2 Definition of * operators (field,constance)
**************************************************************************************/
template<typename T> 
nuc3d::field3d<T>& nuc3d::field3d<T>::operator*= (const T &val)
{
	if (!(field_data==nullptr))
	{
	for(int k=0;k!=nz;++k)
		for(int j=0;j!=ny;++j)
			for(int i=0;i!=nx;++i)
			{
				int idx_loc=getIndex(i,j,k);
				field_data[idx_loc]*=val;
			}
	}
	else
	{
		std::cout<<"ERROR: \n"<<std::endl;
		std::cout<<"Null field for operation '*=' !\n"<<std::endl;		
	}

	return *this;
};

template<typename T> 
nuc3d::field3d<T> nuc3d::field3d<T>::operator* (const T &val)
{
	field3d tmp(*this);
	tmp*=val;

	return tmp;
};

/**************************************************************************************
				2. Definition of / operators (field,field),(field,constance)
**************************************************************************************/
/**************************************************************************************
				2.1 Definition of / operators (field,field)
**************************************************************************************/
template<typename T> 
nuc3d::field3d<T>& nuc3d::field3d<T>::operator/= (const field3d<T> &rhs)
{
	if (nx==rhs.getSizeX()&&
		ny==rhs.getSizeY()&&
		nz==rhs.getSizeZ())
	{
		if (!(field_data==nullptr)&&!(rhs.field_data==nullptr))
		{
			for(int k=0;k!=nz;++k)
				for(int j=0;j!=ny;++j)
					for(int i=0;i!=nx;++i)
					{
						int idx_loc=getIndex(i,j,k);
						field_data[idx_loc]/=rhs.field_data[idx_loc];
					}
		}
		else
		{
			std::cout<<"ERROR: \n"<<std::endl;
			std::cout<<"Null field for operation '/=' !\n"<<std::endl;		
		}
	}
	else
	{
		std::cout<<"ERROR: \n"<<std::endl;
		std::cout<<"Field size does not match for operation '/=' !\n"<<std::endl;		
	}

	return *this;
};

template<typename T> 
nuc3d::field3d<T> nuc3d::field3d<T>::operator/ (const field3d<T> &rhs)
{
	field3d tmp(*this);
	tmp/=rhs;

	return tmp;
};

/**************************************************************************************
				2.2 Definition of / operators (field,constance)
**************************************************************************************/
template<typename T> 
nuc3d::field3d<T>& nuc3d::field3d<T>::operator/= (const T &val)
{
	if (!(field_data==nullptr)&&!(val==0.0))
	{
		for(int k=0;k!=nz;++k)
			for(int j=0;j!=ny;++j)
				for(int i=0;i!=nx;++i)
				{
					int idx_loc=getIndex(i,j,k);
					field_data[idx_loc]/=val;
				}
	}
	else
	{
		std::cout<<"ERROR: \n"<<std::endl;
		if (field_data) 
			std::cout<<"Singular field A for operation '(field A)/=val' !\n"<<std::endl;
		if (val==0.0) 				
			std::cout<<"Singular val for operation '(field A)/=val' !\n"<<std::endl;
	}
	return *this;
};

template<typename T> 
nuc3d::field3d<T> nuc3d::field3d<T>::operator/ (const T &rhs)
{
	field3d tmp(*this);
	tmp/=rhs;

	return tmp;
};

#endif
/**************************************************************************************
									End of this file
**************************************************************************************/