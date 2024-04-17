#pragma once

#include "numerical.hpp"
#include "invparam.hpp"
#include "shared/csr_matrix.hpp"
#include "SWD/SurfaceWave.hpp"

class JSurfGParams: public InverseParamsBase
{
public:
    float noilevel1,noilevel2; // std of gaussion distribution
    float weight1, weight2; // weight for each dataset
    float p; //  relative weight params, see julia et al(2000)
    int rm_grav_avg;  // remove gravity average value 
};

class JSurfGTomo{
public:
    //swd data
    SurfTime surf;

    // gravity data
    csr_matrix gmat;
    fvec obsg,lon_grav,lat_grav; // observed gravity and system

    // inv params
    JSurfGParams param;

    // models
    fmat3 vstrue,vsinit,vsref;
    fvec lon,lat,dep;

private:
    void compute_grav_grad(const fmat3 &vs,fvec &dgsyn,fvec &grad) const;
    
public:
    void read_model(const std::string &modfile,const std::string &modtrue,
                    const std::string &modref);
    void read_data(const std::string &swdfile,const std::string &gravfile);
    void read_gravmat(const std::string &gravfile);
    void read_invparams(const std::string &paramfile);
    void checkerboard();
    void compute_gravity(const fmat3 &vs,fvec &dgsyn) const;
    void inversion(fmat3 &vsf,fvec &dsyn) const;
    void write_syn(const fvec &dysn,const std::string &swdfile,
                    const std::string &gravfile) const;

    // for cg 
    float compute_misfit(const fvec &dsyn) const;
    void forward(const fvec &x,fvec &dsyn) const;
    void compute_grad(const fvec &x,fvec &dsyn,fvec &grad) const;
    void smoothing(fvec &grad) const;
};