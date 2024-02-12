#pragma once 

#include <string>
#include <array>

class InverseParamsBase
{
public:
    int inv_method;
    float minvel,maxvel; // velocity prior information
    int maxiter; // maxiteration
    int ifsyn; // if or not do checkerboard test

    // lsmr
    float smooth,damp; // parameters for lsmr
    int nthreads;

    // nonlinear cg/lbfgs
    int smooth_in_km;
    float sigma_h,sigma_v; // cg smoothing parameters

    // model flags
    int iter_cur,iter_start;

    // line search parameters
    float MAX_REL_STEP;

public:
    void read_file(const std::string &parafile);
};