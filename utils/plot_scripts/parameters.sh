
# path to results directory
RESULT_DIR=../../example/results/

# python scripts
PY_SCRIPTS=./src/

# horizontal slice 
NSLICE_H=1
DEPTH_H=(5 25)
LON0_H=99.7; LON1_H=110.2
LAT0_H=25.7; LAT1_H=35.3

# vertical slice
NSLICE_V=2
MIN_DEPTH=0.
MAX_DEPTH=80.
LON0_V=(105 99.7)
LON1_V=(105 110.2)
LAT0_V=(35.3 30.)
LAT1_V=(25.7 30.)

# options
SYN_TEST=0
HAS_TOPO=0

# grid size
NGRD=256

# dispersion file
DISP_FILE=../../example/surfdataSC.dat

# model used
run_idx="0 10"
