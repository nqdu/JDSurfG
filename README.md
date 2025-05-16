# Joint Inversion of Direct Surface Wave Tomography and Bouguer Gravity
Version 4.2

# Futures
From version 5.0, `Vp` and `Vs` can be independently handled.

## Note
Please note that the format of the input files after `v4.0` are a little different from previous versions!

## New Features
* More accurate group velocity ray-tracing implementation in `fmst` package.
* Add nonlinear conjugate gradient and L-BFGS optimizer for large-scale problems
* High order dispersion data can be included.

## Updates
* Move user manual here [doc](doc/UserManual.md)
* optimized L-BFGS framework, Wolfe conditions are applied.
* Add `clang` support on MacOS.