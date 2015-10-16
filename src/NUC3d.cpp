# include "nuc3d.h"
#include "IOController.h"
#include "field.h"
#include <vector>
#define N (8)
int main ( int argc, char *argv[] )
{
  	nuc3d::IOController myIO;

  	std::cout<<myIO.getMethod("endStep")<<std::endl;
  	std::cout<<myIO.getValue("cfl")<<std::endl;
 	
 	std::fstream myFile;
 	std::vector<nuc3d::field3d<double>> blocks;
 	std::vector<nuc3d::field3d<double>> blocks2(3,nuc3d::field3d<double>(N,N,N,1.0));
 	blocks.push_back(nuc3d::field3d<double>(N,N,N,1.0));
 	blocks.push_back(nuc3d::field3d<double>(N,N,N,2.0));
 	blocks.push_back(nuc3d::field3d<double>(N,N,N,3.0));
	int nx,ny,nz;
	int nx0,ny0,nz0;
	nx=N;
	ny=N;
	nz=N;
	std::stringstream inout;
	std::string filename;
	inout<<"mesh_"<<N<<".dat";
	inout>>filename;
	myFile.open(filename,std::fstream::out);
	//myFile.open(filename,ifstream::binary);
	if(myFile)
	{		
		myFile<<nx<<" "<<ny<<" "<<nz<<std::endl;
		for(int k=0;k<N;k++)
			for(int j=0;j<N;j++)
				for(int i=0;i<N;i++)
				{
					myFile.precision(16);
					for(auto iter=blocks2.begin();iter!=blocks2.end();iter++)
						myFile<<std::fixed<<std::showpoint<<iter->getValue(i,j,k)<<" ";
					myFile<<std::endl;
				}
	}
	else
	{
		std::cout<<"file "<<filename<<" does not exist"<<std::endl;
	}
	myFile.close();

	myFile.open(filename,std::fstream::in);

	if(myFile)
	{
		myFile>>nx0>>ny0>>nz0;
		std::cout<<nx0<<" "<<ny0<<" "<<nz0<<" "<<std::endl;

	}
	std::cout<<blocks.size()<<std::endl;
  return 0;
}

 // nuc3d::MPIComunicator3d_nonblocking<double> myComm(argc,argv);
 //  nuc3d::field3d<double> myField1(N,N,N,myComm.getMyId());
 //  nuc3d::field3d<double> myField2(N,N,N,myComm.getMyId()+1.0);
 //  nuc3d::field3d<double> myField3(N,N,N,myComm.getMyId()+2.0);
 //  std::vector<nuc3d::field3d<double>> Q_Euler(6,nuc3d::field3d<double>(N,N,N));
 //  nuc3d::bufferData<double> myBuffer1;  

 //  std::cout<<"Vector size is "<<Q_Euler.size()<<std::endl;

 //  int myNeibs[2][6]={{-1,1,-1,-1,-1,-1},
 //                     {0,-1,-1,-1,-1,-1}};
 //  int myNeibFaces[2][6]={{-1,0,-1,-1,-1,-1},
 //                    {1,-1,-1,-1,-1,-1}};

 //  myComm.setTopo(myNeibs[myComm.getMyId()],myNeibFaces[myComm.getMyId()]);
 //  myBuffer1.allocatebuffer(N,N,N,3);

 //  myBuffer1.setBufferSend(myField1);
 //  myComm.bufferSendRecv(myBuffer1);
 //  myComm.waitAllSendRecv();

 //  if(myComm.getMyId()==0)
 //    std::cout<<"buffer data for id "<<myComm.getMyId()<<" is "
 //            <<myBuffer1.BufferRecv[1][0]<<"\n"
 //            <<myField1.getValue(0,0,0)<<"\n"
 //            <<std::endl;


 //  myBuffer1.setBufferSend(myField2);
 //  myComm.bufferSendRecv(myBuffer1);
 //  myComm.waitAllSendRecv();

 //  if(myComm.getMyId()==0)
 //    std::cout<<"buffer data for id "<<myComm.getMyId()<<" is "
 //            <<myBuffer1.BufferRecv[1][0]<<"\n"
 //            <<myField2.getValue(0,0,0)<<"\n"
 //            <<std::endl;

 //  myBuffer1.setBufferSend(myField3);
 //  myComm.bufferSendRecv(myBuffer1);
 //  myComm.waitAllSendRecv();

 //  if(myComm.getMyId()==0)
 //    std::cout<<"buffer data for id "<<myComm.getMyId()<<" is "
 //            <<myBuffer1.BufferRecv[1][0]<<"\n"
 //            <<myField3.getValue(0,0,0)<<"\n"
 //            <<std::endl;

 //  myBuffer1.deallocatebuffer();
  //MPI_Barrier(MPI_COMM_WORLD);