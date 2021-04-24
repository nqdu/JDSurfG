# compiler and options
include ../../include/Makefile
CC = g++
FC = gfortran
FFLAGS  = -ffixed-line-length-none -ffloat-store\
           -w -m64 -mcmodel=medium
includes = -I../../include $(path_to_eigen)

# c scripts
src = SurfaceWave.cpp 

# fortran src
f90mp=fmm.f90 fmst.f90 slegn96.f90 sregn96.f90
f77src = surfdisp96.f

# objective
cppobj = $(src:%.cpp=%.o) surfdisp.o 
f90obj= $(f90mp:%.f90=%.o) 
f77obj=$(f77src:%.f=%.o) 

main:$(cppobj) $(f90obj) $(f77obj)

%.o:%.f90
	$(FC) -g -c $< -o $@  $(opt) -fopenmp

%.o:%.cpp
	$(CC) -g -c $(cppstd) $(includes)  $< -o $@ $(opt) -fopenmp 

surfdisp.o:surfdisp.cpp
	$(CC) -g -c $(cppstd) $(includes)  $< -o $@ $(opt)

%.o:%.f
	$(FC) $(FFLAGS) -g -rdynamic -c $< -o $@ $(opt)	


clean:
	rm  *.o *.mod