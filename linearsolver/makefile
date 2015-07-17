objects = solver.o basic_operation.o convert.o read.o serial_solve.o xvalue.o
CXX = icc
CXXFLAGS = -openmp 
solve : $(objects)
	icc -o solve $(objects) -openmp
solver.o : basic_operation.h IO_functions.h parameters.h

basic_opration.o : IO_functions.h parameters.h

convert.o : parameters.h

read.o : IO_functions.h parameters.h

serial_solve.o : basic_operation.h IO_functions.h parameters.h

xvalue.o : IO_functions.h basic_operation.h

.PHONY : clean
clean :
	rm solve $(objects)
