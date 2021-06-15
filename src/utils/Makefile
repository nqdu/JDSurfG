# compiler and options
include ../../include/Makefile
includes = -I../../include

# c scripts
src = gaussian.cpp spherical.cpp

# objective
cppobj = $(src:%.cpp=%.o) gauleg.o bkg_model.o IOFunction.o
f90obj= lsmrModule_csr.o 

main: $(f90obj)  $(cppobj)

lsmrModule_csr.o : lsmrModule_csr.f90
	$(FC) -g -c $< -o $@  $(opt) -fopenmp

%.o:%.cpp
	$(CC) $(CXXFLAGS) -g -c $(includes)  $< -o $@

gauleg.o:gauleg.cpp 
	$(CC) $(CXXFLAGS) -g -c $(includes) $(GSL_INC)  $< -o $@ 

IOFunction.o:IOFunction.cpp 
	$(CC) $(CXXFLAGS) -g -c $(cppstd) $(includes) $(EIGEN_INC)   $< -o $@ 
bkg_model.o:bkg_model.cpp
	$(CC) $(CXXFLAGS) -g -c $(cppstd) $(includes) $(EIGEN_INC)   $< -o $@ 

clean:
	rm  *.o *.mod