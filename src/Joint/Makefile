out = JDsurfG

# compiler and options
opt=-O3
cc = g++
FC = gfortran
FFLAGS  = -ffixed-line-length-none -ffloat-store\
           -w  -fbounds-check -m64 -mcmodel=medium
includes = -I../../include

# c scripts
src = $(shell find *.cpp) ../defmod.cpp ../spherical.cpp ../mod3d.cpp ../surftime.cpp

# fortran src
f90src = ../gaussian.f90 ../lsmrModule.f90 ../fmm.f90  ../delsph.f90 ../empirical.f90
f77src = ../surfdisp96.f

# objective
cppobj = $(src:%.cpp=%.o) 
fobj=  $(f90src:%.f90=%.o) $(f77src:%.f=%.o) ../calsurf.o

$(out): $(fobj) $(cppobj)
	$(cc)  -o ../../bin/$(out) $(fobj) $(cppobj) -lgfortran -fopenmp

%.o: %.f90
	$(FC) -g -c $< -o $@  $(opt)

../calsurf.o:../calsurf.f90
	$(FC) -g -c $< -o $@  $(opt) -fopenmp

%.o: %.f
	$(FC) $(FFLAGS) -g -rdynamic -c $< -o $@ $(opt)	

%.o: %.cpp $(deps)
	$(cc) -g -c -std=c++11 $(includes) $< -o $@ $(opt)

clean:
	rm  *.o 
	