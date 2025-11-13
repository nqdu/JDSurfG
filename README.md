# Joint Inversion of Direct Surface Wave Tomography and Bouguer Gravity
Version 4.3

# Future work
From version 5.0, `Vp` and `Vs` can be independently handled.

## Note
From version 4.3, we migrate from `C++11` to `C++14`. 

## New Features
* Add topographic correction
* More accurate group velocity ray-tracing implementation in `fmst` package.
* Add nonlinear conjugate gradient and L-BFGS optimizer for large-scale problems
* High order dispersion data can be included.

## Updates
* Move user manual here [doc](doc/UserManual.md)
* optimized L-BFGS framework, Wolfe conditions are applied.
* Add `clang` support on MacOS.