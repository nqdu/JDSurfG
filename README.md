# Joint Inversion of Direct Surface Wave Tomography and Bouguer Gravity
Version 4.3

# Future work
From version 5.0, `Vp` and `Vs` can be independently handled.

## New Features
* Add topographic correction
* Accurate group velocity ray-tracing implementation in `fmst` package.
* Nonlinear conjugate gradient and L-BFGS optimizer for large-scale problems
* High order dispersion data can be included.

## Updates
* Move user manual here [doc](doc/UserManual.md)
* optimized L-BFGS framework, Wolfe conditions are applied.
* Add `clang` support on MacOS.