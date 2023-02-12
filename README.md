# Joint Inversion of Direct Surface Wave Tomography and Bouguer Gravity
Version 3.2. Update from version 3.0. Version 3.0+ is about **5-10 times faster** than previous versions!

## Note
Please note that the format of the input files are a little different from previous versions!

## New Features
* Parallel linear system solvers in `src/utils/lsmrModule_csr.f90`.
* Rewrite `slegn96` and `sregn96` in [Computer Programs in Seismology](http://www.eas.slu.edu/eqc/eqccps.html). Now analytical derivatives for Love and Rayleigh wave dispersion could be utilized.

## Updates
* Update parameter files, docs and references.
* Add I/O functions.
* Add intel compiler support 

## Bug Fix
* Fix bugs in `utils/gravity_forward.f90`
* Fix bugs in the computing gravity matrix.

## TODO List (from easy to hard)  


~~1. A more concise representation of station pairs with C++ class.~~

~~2. Parallelization of the computation of surface wave travetime and 2-D Frechet kernel.~~

~~3. Inhomegeneous grid in vertical direction.~~

~~4. Analytical sensitvity kernels for surface wave dispersions~~

~~5. Water layers could be involved in inversion.~~ 

6. Different empirical relations for different region

7. Utilization of higher mode surface wave dispersion.

