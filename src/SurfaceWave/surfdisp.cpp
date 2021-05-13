#include"surfdisp.hpp"
#include<math.h>

/**
 * convert phase/group velocity of flat earth to that of spherical earth. 
 * Refer to Schwab, F. A., and L. Knopoff (1972). Fast surface wave and free
 * mode computations, in  Methods in Computational Physics, Volume 11,
 * Love Wave Equations  44, 45 , 41 pp 112-113. 
 * Rayleigh Wave Equations 102, 108, 109 pp 142, 144.
 * @param t period
 * @param c phase/group velo
 * @param wavetp wavetype one of [Rc,Rg,Lc,Lg]
 * 
 * @return cnew phase/group velocity of spherical earth
 */
double _flat2sphere(double t,double c,std::string wavetp)
{
    double ar = 6371.0; // earth radius
    double tm;
    double omega = 2.0 * M_PI / t; // angular frequency

    // check input parameters
    bool flag = wavetp == "Rc" || wavetp == "Rg" || 
                wavetp == "Lc" || wavetp == "Lg";
    if(flag == false){
        std::cout << "cnew = _flat2sphere(double t,double c,std::string wavetp)\n";
        std::cout << "parameter wavetp should be in one of [Rc,Rg,Lc,Lg]\n ";
        exit(0);
    }

    if(wavetp[0] == 'L'){ //love wave
        tm = 1. + pow(1.5 * c / (ar * omega),2); 
    }
    else{ // Rayleigh wave
        tm = 1. + pow(0.5 * c / (ar * omega),2);
    }
    tm = sqrt(tm);

    // convert to spherical velocity
    double cnew;
    if(wavetp[1] == 'c'){ // phase velocity
        cnew = c / tm;
    }
    else{ // group velocity
        cnew = c * tm;
    }

    return cnew;
}

/**
 * calculates the dispersion values for any layered model, any frequency, and any mode.
 * @param nlayer no. of layers
 * @param thkm,vpm,vsm,rhom model, shape(nlayer)
 * @param kmax no. of periods used
 * @param t,cp period and phase velocity, shape(kmax)
 * @param sphere true for spherical earth, false for flat earth 
 * @param wavetype one of [Rc,Rg,Lc,Lg]
 * @param mode i-th mode of surface wave, 0 fundamental, 1 first higher, ....
 * @param keep_flat keep flat earth phase/group velocity or convert it to spherical
 */
void surfdisp(float *thk,float *vp,float *vs,float *rho,
            int nlayer,double *t,double *cg,int kmax,const std::string &wavetype,
            int mode,bool sphere,bool keep_flat)
{
    int iwave,igr,ifsph=0;
    if(wavetype=="Rc"){
        iwave = 2;
        igr = 0;
    }
    else if(wavetype == "Rg"){
        iwave = 2;
        igr = 1;
    }
    else if(wavetype=="Lc"){
        iwave = 1;
        igr = 0;
    }
    else if(wavetype=="Lg"){
        iwave = 1;
        igr = 1;
    }
    else{
        std::cout <<"wavetype should be one of [Rc,Rg,Lc,Lg]"<<std::endl;
        exit(0);
    }

    if(sphere == true) ifsph = 1;
    surfdisp96_(thk,vp,vs,rho,nlayer,ifsph,iwave,mode+1,igr,kmax,t,cg);

    if(sphere == true && keep_flat == false){
        for(int i=0;i<kmax;i++){
            cg[i] = _flat2sphere(t[i],cg[i],wavetype);
        }
    }
}

/** Love group velocity with analytical method  
 * @param nlayer no. of layers
 * @param thk,vs,rho 1-D model
 * @param kmax no. of periods
 * @param t cg periods and group velocity
 * @param mode i-th mode of surface wave, 1 fundamental, 2 first higher, ....
 * @param sphere true for spherical earth, false for flat earth 
 */
void _LoveGroup(float *thk,float *vs,float *rho,
              int nlayer,double *t,double *cg,int kmax,
              int mode,bool sphere)
{
    // check if earth model is spherical
    int iflsph = 0; if(sphere == true) iflsph = 1;

    // temporay arrays
    float vp[kmax];
    double cp[kmax],uu[nlayer],tt[nlayer],dcdh[nlayer],
            dcdr[nlayer],dcdb[nlayer];

    // compute phase velocity
    for(int i=0;i<kmax;i++) vp[i] = 1.732 * vs[i];
    surfdisp(thk,vp,vs,rho,nlayer,t,cp,kmax,"Lc",mode,sphere,true);

    // convert phase to group velocity
    for(int i=0;i<kmax;i++){
        slegn96_(thk,vs,rho,nlayer,t+i,cp+i,cg+i,uu,tt,dcdb,dcdh,dcdr,iflsph);
    }
}

/** Rayleigh group velocity with analytical method  
 * @param nlayer no. of layers
 * @param thk,vp,vs,rho 1-D model
 * @param kmax no. of periods
 * @param t cg periods and group velocity
 * @param mode i-th mode of surface wave, 1 fundamental, 2 first higher, ....
 * @param sphere true for spherical earth, false for flat earth 
 */
void _RayleighGroup(float *thk,float *vp,float *vs,float *rho,
              int nlayer,double *t,double *cg,int kmax,
              int mode,bool sphere)
{
    // check if earth model is spherical
    int iflsph = 0; if(sphere == true) iflsph = 1;

    // temporay arrays
    double cp[kmax],ur[nlayer],uz[nlayer],tr[nlayer],tz[nlayer],
            dcdh[nlayer],dcda[nlayer],dcdr[nlayer],dcdb[nlayer];

    // compute phase velocity
    surfdisp(thk,vp,vs,rho,nlayer,t,cp,kmax,"Rc",mode,sphere,true);

    // convert phase to group velocity
    for(int i=0;i<kmax;i++){
        sregn96_(thk,vp,vs,rho,nlayer,t+i,cp+i,cg+i,
                ur,uz,tr,tz,dcda,dcdb,dcdh,dcdr,iflsph);
    }
}