out = DSurfTomo

# compiler and options
include ../../include/Makefile
cc = g++
FC = gfortran
FFLAGS  = -ffixed-line-length-none -ffloat-store\
           -w -m64 -mcmodel=medium
includes = -I../../include $(path_to_eigen)
surfdir= ../SurfaceWave

# c scripts
src = surftomo.cpp main.cpp $(surfdir)/SurfaceWave.cpp \
		../utils/gaussian.cpp ../utils/bkg_model.cpp

# fortran src
f90mp=$(surfdir)/fmm.f90 $(surfdir)/fmst.f90 
f90src =  ../utils/delsph.f90 ../utils/empirical.f90 ../utils/lsmrModule_csr.f90
f77src = $(surfdir)/surfdisp96.f

# objective
cppobj = $(src:%.cpp=%.o)
fobj= $(f90src:%.f90=%.o) $(f90mp:%.f90=%.o) $(f77src:%.f=%.o) 

$(out): $(fobj) $(cppobj)
	$(cc)  -o ../../bin/$(out) $(fobj) $(cppobj) -lgfortran -fopenmp

%.o: %.cpp
	$(cc) -g -c $(cppstd) $(includes) $< -o $@ $(opt)

clean:
	rm  *.o 
