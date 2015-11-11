# include "NUC3d.h"

int main ( int argc, char *argv[] )
{
    if(MPI_Init( &argc, &argv ))
    {
        std::cout<<"ERROR: MPI Initialization failed!"<<std::endl;
        MPI_Abort ( MPI_COMM_WORLD,99 );
    }
    
    nuc3d::singleBlock myBlock;    
    return 0;
}
