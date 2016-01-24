 cc=mpicxx
 srcs=block.cpp\
boundaryConditions.cpp\
bufferData.cpp\
centraldifference6th.cpp\
centraldifference8th.cpp\
euler3d.cpp\
NaiverStokes3d.cpp\
field.cpp\
fieldOperator.cpp\
IOcontroller.cpp\
MPICommunicatorBase.cpp\
MPICommunicator.cpp\
NUC3d.cpp\
PDEData3d.cpp\
physicsModel.cpp\
singleBlock.cpp\
TVD-RK3.cpp\
	weno5js.cpp\
	weno7js.cpp\
	hccs.cpp\
 	weno5z.cpp\
	weno3js.cpp\
	crweno5.cpp\
	upwind1st.cpp\
	centraldifference2nd.cpp\
	boundaryScheme.cpp\
	gradvector.cpp\
	reconstructionBoundaryScheme.cpp\
	postproc.cpp
OBJS=$(srcs:.cpp=.o)

%.o:%.cpp
	$(cc) -std=c++11 -c $<
default:$(OBJS)
	$(cc) -O2 -o a.out $(OBJS)
install:
	cp ./a.out ../../bin
clean:
	rm -f *.out *.o
