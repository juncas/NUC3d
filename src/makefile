 cc=mpicxx
 srcs=block.cpp\
boundaryConditions.cpp\
bufferData.cpp\
centraldifference6th.cpp\
euler3d.cpp\
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
 	weno5z.cpp\
	upwind1st.cpp
OBJS=$(srcs:.cpp=.o)

%.o:%.cpp
	$(cc) -std=c++11 -c $<
default:$(OBJS)
	$(cc) -O2 -o a.out $(OBJS)
clean:
	rm -f *.out *.o
