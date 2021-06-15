src=$(shell find *.f90)
obj = $(src:%.f90=%.o)
out = ../bin/syngrav
include ../include/Makefile

$(out):  $(obj)
	$(FC) -o $(out) $(obj) 

%.o: %.f90
	$(FC) $(FFLAGS) -g -c $< -o $@ 

clean:
	rm *.o *.mod
