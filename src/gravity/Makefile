out = mkmat

# compiler and options
opt=-O3
cc = g++
includes = -I../../include
gsl_lib = -lgsl -lgslcblas

# c scripts
src =../defmod.cpp  gauleg.cpp  main.cpp ../spherical.cpp

# objective
cppobj = $(src:%.cpp=%.o) gravmat.o ../empirical.o 

$(out):  $(cppobj)
	$(cc)  -o ../../bin/$(out) $(cppobj) -lm -fopenmp $(gsl_lib)

gravmat.o:gravmat.cpp
	$(cc) -g -c -std=c++11 $(includes) $< -o $@ $(opt) -fopenmp
%.o: %.cpp
	$(cc) -g -c -std=c++11 $(includes) $< -o $@ $(opt)

%.o: %.f90
	gfortran -g -c $< -o $@ $(opt)
clean:
	rm  *.o