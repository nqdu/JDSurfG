# compiler and FLAGS
# GNU
CC = g++ 
FC = gfortran 
CXXFLAGS = -std=c++11 -O3 -march=native
FFLAGS = -O3 -march=native -ffixed-line-length-none
LINK_CF = -lgfortran

#intel
# CC = icc 
# FC = ifort
# CXXFLAGS = -std=c++11 -O3 -xHost -extend-source
# FFLAGS = -O3 -xHost 
# LINK_CF = -lifcore

# library path
#EIGEN_INC = -I/home/dunanqiao/software/eigen-3.3.9
EIGEN_INC = -I/mnt/c/linuxapp/eigen-3.4.0
GSL_LIB = $(shell gsl-config --libs)
GSL_INC = $(shell gsl-config --cflags)

