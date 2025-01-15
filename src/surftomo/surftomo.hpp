#ifndef JDSURFG_SURFTOMO_SURFTOMO_H_
#define JDSURFG_SURFTOMO_SURFTOMO_H_


#include "numerical.hpp"
#include "invparam.hpp"
#include "SWD/SurfaceWave.hpp"

class DSurfTomoParams: public InverseParamsBase
{
public:
    float noiselevel; // std of gaussion distribution
};


class DSurfTomo{
public:
    SurfTime surf;
    DSurfTomoParams param;
    fmat3 vstrue,vsinit;
    fvec lon,lat,dep;
    
public:
    void read_model(const std::string &modfile,const std::string &modtrue) ;
    void read_data(const std::string &datafile);
    void read_invparams(const std::string &paramfile);
    void checkerboard();
    void inversion(fmat3 &vsf,fvec &dsyn);

    // for cg 
    float compute_misfit(const fvec &dsyn) const;
    void forward(const fvec &x,fvec &dsyn) const;
    void compute_grad(const fvec &x,fvec &dsyn,fvec &grad) const;
    void smoothing(fvec &grad) const;
};

#endif // end JDSURFG_SURFTOMO_SURFTOMO_H_
