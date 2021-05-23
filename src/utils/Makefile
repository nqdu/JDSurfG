# compiler and options
include ../../include/Makefile
cc = g++
FC = gfortran
includes = -I../../include

# c scripts
src = gaussian.cpp spherical.cpp

# objective
cppobj = $(src:%.cpp=%.o) gauleg.o bkg_model.o
f90obj= lsmrModule_csr.o 

main: $(f90obj)  $(cppobj)

lsmrModule_csr.o : lsmrModule_csr.f90
	$(FC) -g -c $< -o $@  $(opt) -fopenmp

%.o:%.cpp
	$(cc) -g -c $(cppstd) $(includes)  $< -o $@ $(opt)

gauleg.o:gauleg.cpp 
	$(cc) -g -c $(cppstd) $(includes) $(path_to_gsl_include)  $< -o $@ $(opt)

bkg_model.o:bkg_model.cpp
	$(cc) -g -c $(cppstd) $(includes) $(path_to_eigen)  $< -o $@ $(opt)

clean:
	rm  *.o *.mod