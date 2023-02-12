# compiler and options
include ../../include/Makefile
includes = -I../../include $(EIGEN_INC)
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
	$(FC) $(FFLAGS) -g -c $< -o $@  -fopenmp

%.o:%.cpp
	$(CC) -g -c $(CXXFLAGS) $(includes)  $< -o $@  -fopenmp 

surfdisp.o:surfdisp.cpp
	$(CC) -g -c $(CXXFLAGS) $(includes)  $< -o $@ 

%.o:%.f
	$(FC)  $(FFLAGS)  -g  -c $< -o $@ -fopenmp

clean:
	rm  *.o *.mod