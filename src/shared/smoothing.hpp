#pragma once 
#include "numerical.hpp"
void smooth_cart_pde(float* grad,int nx,int ny,int nz,float sigma_h,float sigma_v);
void smooth_sph_pde(float* gradin,int nx,int ny,int nz,
                float dx,float dy,float dz,
                float lat0,float lon0,float z0,
                float sigma_h,float sigma_v);