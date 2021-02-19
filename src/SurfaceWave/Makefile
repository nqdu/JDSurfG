# compiler and options
include ../../include/Makefile
cc = g++
FC = gfortran
FFLAGS  = -ffixed-line-length-none -ffloat-store\
           -w -m64 -mcmodel=medium
includes = -I../../include $(path_to_eigen)

# c scripts
src = SurfaceWave.cpp 

# fortran src
f90mp=fmm.f90 fmst.f90 
f77src = surfdisp96.f

# objective
cppobj = $(src:%.cpp=%.o)
f90obj= $(f90mp:%.f90=%.o) 
f77obj=$(f77src:%.f=%.o) 

main:$(cppobj) $(f90obj) $(f77obj)

%.o:%.f90
	$(FC) -g -c $< -o $@  $(opt) -fopenmp

%.o:%.cpp
	$(cc) -g -c $(cppstd) $(includes)  $< -o $@ $(opt) -fopenmp 

%.o:%.f
	$(FC) $(FFLAGS) -g -rdynamic -c $< -o $@ $(opt)	


clean:
	rm  *.o *.mod