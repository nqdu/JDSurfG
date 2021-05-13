#pragma once
#include<iostream>
extern "C"{

/**
 * compute phase velocity sensitivity kernel, and group velocity 
 * for Rayleigh Wave,with layer-based model. 
 * note that if in spherical model, input cp should be cp_flat
 * @param nlayer no. of layers
 * @param thk,vp,vs,rhom model, shape(nlayer)
 * @param t,cp,cg period and phase/group velocity
 * @param iflsph =0 flat earth   =1 spherical earth 
 * @param dispu/w stressu/w eigen function for u/w direction shape(nlayer)
 * @param dc2da(b,h,r) sensitivity kernel for vp,vs,thick and rho shape(nlayer)
 *
 */
void sregn96_(float *thk,float *vp,float *vs,float *rhom,
             int nlayer,double *t,double *cp,double *cg,double *dispu,
             double *dispw,double *stressu,double *stressw,
             double *dc2da,double *dc2db,double *dc2dh,
             double *dc2dr,int iflsph);

/**
 * compute phase velocity sensitivity kernel, and group velocity  
 * for Love Wave,with layer-based model
 * note that if in spherical model, input cp should be cp_flat
 * @param nlayer no. of layers
 * @param thk,vs,rhom model, shape(nlayer)
 * @param t,cp,cg period and phase/group velocity
 * @param iflsph =0 flat earth   =1 spherical earth 
 * @param disp,stress eigen function for u/w direction shape(nlayer)
 * @param dc2db(h,r) phase v sensitivity kernel for vs,thick and rho shape(nlayer)
 */
void slegn96_(float *thk,float *vs,float *rhom,int nlayer,
                  double *t,double *cp,double *cg,double *disp,
                  double *stress,double *dc2db,double *dc2dh,
                  double *dc2dr,int iflsph);

/**
 * calculates the dispersion values for any layered model, any frequency, and any mode.
 * @param nlayer no. of layers
 * @param thkm,vpm,vsm,rhom model, shape(nlayer)
 * @param kmax no. of periods used
 * @param t,cp period and phase velocity, shape(kmax)
 * @param iflsph 0 for flat earth and 1 for spherical earth 
 * @param iwave 1 for Love and 2 for Rayleigh 
 * @param mode i-th mode of surface wave, 1 fundamental, 2 first higher, ....
 * @param igr 0 phase velocity, > 0 group velocity
 */
void surfdisp96_(float *thkm,float *vpm,float *vsm,float *rhom,
                int nlayer,int iflsph,int iwave,int mode,
                int igr,int kmax,double *t,double *cg);

/**
 * compute phase/group velocity sensitivity kernel, and group velocity  
 * for Love Wave,with layer-based model
 * note that if in spherical model, input cp should be cp_flat
 * @param nlayer no. of layers
 * @param thk,vs,rhom model, shape(nlayer)
 * @param t,cp,cg period and phase/group velocity
 * @param t1,t2,cp1,cp2 slightly different period and phasev from t
 * @param iflsph =0 flat earth   =1 spherical earth 
 * @param disp,stress eigen function for u/w direction shape(nlayer)
 * @param dc2db(h,r) phase velocity sensitivity kernel for vs,thick and rho shape(nlayer)
 * @param du2db(h,r) group velocity sensitivity kernel for vs,thick and rho shape(nlayer)
 */
void slegnpu_(float *thk,float *vs,float *rhom,int nlayer,
                double *t,double *cp,double *cg,double *disp,
                double *stress,double *t1,double *cp1,double *t2,double *cp2,
                double *dc2db,double *dc2dh,double *dc2dr,
                double *du2db,double *du2dh,double *du2dr,
                int iflsph);

/**
 * compute phase velocity sensitivity kernel, and group velocity 
 * for Rayleigh Wave,with layer-based model. 
 * note that if in spherical model, input cp should be cp_flat
 * @param nlayer no. of layers
 * @param thk,vp,vs,rhom model, shape(nlayer)
 * @param t,cp,cg period and phase/group velocity
 * @param t1,t2,cp1,cp2 slightly different period and phasev from t
 * @param iflsph =0 flat earth   =1 spherical earth 
 * @param dispu/w stressu/w eigen function for u/w direction shape(nlayer)
 * @param dc2da(b,h,r) phase velocity sensitivity kernel for vp,vs,thick and rho shape(nlayer)
 * @param du2da(b,h,r) phase velocity sensitivity kernel for vp,vs,thick and rho shape(nlayer)
 */
void sregnpu_(float *thk,float *vp,float *vs,float *rhom,
             int nlayer,double *t,double *cp,double *cg,double *dispu,
             double *dispw,double *stressu,double *stressw,
             double *t1,double *cp1,double *t2,double *cp2,
             double *dc2da,double *dc2db,double *dc2dh,double *dc2dr,
             double *du2da,double *du2db,double *du2dh,double *du2dr,
             int iflsph);
}

// surface wave dispersion c++ wrapper
void surfdisp(float *thk,float *vp,float *vs,float *rho,
            int nlayer,double *t,double *cg,int kmax,const std::string &wavetype,
            int mode=0,bool sphere=false,bool keep_flat=true);


// Love group velocity with analytical method
void _LoveGroup(float *thk,float *vs,float *rho,
              int nlayer,double *t,double *cg,int kmax,
              int mode=0,bool sphere=false);

// Rayleigh group velocity with analytical method
void _RayleighGroup(float *thk,float *vp,float *vs,float *rho,
              int nlayer,double *t,double *cg,int kmax,
              int mode=0,bool sphere=false);             