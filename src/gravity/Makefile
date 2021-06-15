out = mkmat

# compiler and options
include ../../include/Makefile

includes = -I../../include $(EIGEN_INC)
utils=../utils

# c scripts
src = $(shell find *.cpp) $(utils)/spherical.cpp  $(utils)/gauleg.cpp $(utils)/IOFunction.cpp
# objective
cppobj = $(src:%.cpp=%.o)

$(out):  $(cppobj)
	$(CC)  -o ../../bin/$(out) $(cppobj)  -fopenmp $(GSL_LIB)

gravmat.o:gravmat.cpp
	$(CC) -g -c $(CXXFLAGS) $(includes)  $< -o $@ -fopenmp
%.o: %.cpp
	$(CC) -g -c $(CXXFLAGS) $(includes) $< -o $@

%.o: %.f90
	$(FC) -g -c $(FFLAGS) $< -o $@
clean:
	rm  *.o 
