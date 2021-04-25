# Joint Inversion of Direct Surface Wave Tomography and Bouguer Gravity
Version 3.0. Important version! This version is about **5-10 times faster** than previous versions!

## New Features
* Rewrite `slegn96` and `sregn96` in [Computer Programs in Seismology](http://www.eas.slu.edu/eqc/eqccps.html). Now analytical derivatives for Love and Rayleigh wave dispersion could be utilized.
* Water layer could be involved in joint inversion.

## Bug Fix
* Fix the compiled name of function `surfdisp96`  to `surfdisp96_` in `src/SurfaceWave/surfdisp96.f`.
* Fix bugs in parallel reading.
## TODO List (from easy to hard)  


~~1. A more concise representation of station pairs with C++ class.~~

~~2. Parallelization of the computation of surface wave travetime and 2-D Frechet kernel.~~

~~3. Inhomegeneous grid in vertical direction.~~

~~4. Analytical sensitvity kernels for surface wave dispersions~~

~~5. Water layers could be involved in inversion.~~ 

6. Different empirical relations for different region

7. Utilization of higher mode surface wave dispersion.

