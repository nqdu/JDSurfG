out = mkmat

# compiler and options
include ../../include/Makefile

cc = g++
includes = -I../../include $(path_to_eigen) $(path_to_gsl_include)
gsl_lib = $(path_to_gsl_lib)
utils=../utils

# c scripts
src = $(shell find *.cpp) $(utils)/spherical.cpp  $(utils)/gauleg.cpp
# objective
cppobj = $(src:%.cpp=%.o) $(utils)/empirical.o

$(out):  $(cppobj)
	$(cc)  -o ../../bin/$(out) $(cppobj)  -fopenmp $(gsl_lib)

gravmat.o:gravmat.cpp
	$(cc) -g -c $(cppstd) $(includes)  $< -o $@ $(opt) -fopenmp
%.o: %.cpp
	$(cc) -g -c $(cppstd) $(includes) $< -o $@ $(opt)

%.o: %.f90
	gfortran -g -c $< -o $@ $(opt)
clean:
	rm  *.o