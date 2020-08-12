src=$(shell find *.f90)
obj = $(src:%.f90=%.o)
out = ../bin/syngrav

$(out):  $(obj)
	gfortran -o $(out) $(obj) 

%.o: %.f90
	gfortran -g -c $< -o $@ $(opt)

clean:
	rm *.o *.mod