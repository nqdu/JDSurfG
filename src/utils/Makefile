# compiler and options
include ../../include/Makefile
cc = g++
FC = gfortran
includes = -I../../include

# c scripts
src = gaussian.cpp spherical.cpp

# fortran src
f90mp=empirical.f90 delsph.f90

# objective
cppobj = $(src:%.cpp=%.o) gauleg.o bkg_model.o
f90obj= $(f90mp:%.f90=%.o) lsmrModule_csr.o 

main: $(f90obj)  $(cppobj)

%.o:%.f90
	$(FC) -g -c $< -o $@  $(opt)

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