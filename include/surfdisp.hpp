#pragma once
#include<string>
#include<iostream>
extern "C"{
    void surfdisp96(float *thk,float *vp,float *vs,
                    float *rho,int nlayer,int iflsph,int iwave,
                    int mode,int igr,int kmax,double *t,double *cg);
    void synthetic(int nx,int ny,double *velf,float goxdf,float gozdf,
                float dvxdf,float dvzdf,float srcx,float srcz,
                float *rcx,float *rcz,int nr,float *ttime);
    void CalFrechet(int nx,int ny,int nz,double *velf,float goxdf,float gozdf,
                    float dvxdf,float dvzdf,float srcx,float srcz,float *rcx,
                    float *rcz,int nr,float *ttime,double *kernel,float *frechet);
}

void surfdisp(float *thk,float *vp,float *vs,float *rho,
            int nlayer,double *t,double *cg,int kmax,std::string wavetype,
            int mode=0,bool flatearth=false)
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

    if(!flatearth) ifsph = 1;

    surfdisp96(thk,vp,vs,rho,nlayer,ifsph,iwave,mode+1,igr,kmax,t,cg);
}