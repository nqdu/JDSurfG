out = JointSG

# compiler and options
include ../../include/Makefile
includes = -I../../include $(EIGEN_INC) $(GSL_INC)
surfdir= ../SurfaceWave
utils = ../utils

# c++ scripts
src =  main.cpp $(surfdir)/SurfaceWave.cpp $(utils)/gaussian.cpp \
	   $(utils)/spherical.cpp $(utils)/gauleg.cpp ../gravity/gravmat.cpp \
	   $(utils)/bkg_model.cpp $(surfdir)/surfdisp.cpp \
	   $(utils)/IOFunction.cpp


# fortran src
f90mp=$(surfdir)/fmm.f90 $(surfdir)/fmst.f90  \
		$(surfdir)/slegn96.f90  $(surfdir)/sregn96.f90 
f90src =  ../utils/lsmrModule_csr.f90
f77src = $(surfdir)/surfdisp96.f

# objective
cppobj = $(src:%.cpp=%.o) JointTomo.o
fobj= $(f90src:%.f90=%.o) $(f90mp:%.f90=%.o) $(f77src:%.f=%.o) 

$(out): $(fobj) $(cppobj)
	$(CC) -o ../../bin/$(out) $(fobj) $(cppobj) -lgfortran -fopenmp $(GSL_LIB)

%.o: %.cpp
	$(CC) $(CXXFLAGS) -g -c $(cppstd) $(includes) $< -o $@

JointTomo.o:JointTomo.cpp
	$(CC) $(CXXFLAGS) -g -c $(cppstd) $(includes) $< -o $@ -fopenmp
clean:
	rm  *.o
