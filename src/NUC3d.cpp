# include "nuc3d.h"

int main ( int argc, char *argv[] )
{
    if(MPI_Init( &argc, &argv ))
    {
        std::cout<<"ERROR: MPI Initialization failed!"<<std::endl;
        MPI_Abort ( MPI_COMM_WORLD,99 );
    }
    
    nuc3d::singleBlock myBlock;
    
    myBlock.initialBlock();
    myBlock.solvePDE();
    myBlock.postprocess();
    myBlock.output();
    
    MPI_Finalize();
    return 0;
}
