#include <iostream>
#include "SWD/SurfaceWave.hpp"
#include <cmath>

int main (){
    int nx = 33 ,ny = 36;
    float dx = 0.3,dy = 0.3,lat0 = 35.3,lon0 = 99.7;
    fmst_init(nx,ny,lat0,lon0,dx,dy);

    // velocity
    dmat2 vel(nx,ny);
    fmat2 fdm(nx,ny);
    vel.setConstant(3.);

    // source and receiver
    const float deg2rad = M_PI/ 180.;
    float srcx = 30,srcz = 102;
    const int nr = 1;
    float rcx[nr] = {30}, rcz[nr] = {109.};
    float t[nr];

    // print source/receiver
    FILE *fps = fopen("sr.dat","w");
    fprintf(fps,"%f %f\n",srcz,srcx);
    for(int ir = 0; ir < nr; ir ++) {
        fprintf(fps,"%f %f\n",rcz[ir],rcx[ir]);
    }
    fclose(fps);
    
    // convert source/receiver to colat/lon in rad
    srcx = (90 - srcx) * deg2rad;
    srcz = srcz * deg2rad;
    for(int ir = 0; ir < nr; ir ++) {
        rcx[ir] = (90 - rcx[ir]) * deg2rad;
        rcz[ir] = rcz[ir] * deg2rad;
    }

    // phase velocity fdm
    fmst_run(vel.data(),srcx,srcz,rcx,rcz,nr,t);
    fmst_reset(vel.data());
    for(int ir = 0; ir < nr; ir ++) {
        fmst_raypath(srcx,srcz,rcx[ir],rcz[ir],fdm.data());
    }

    // write output 
    FILE *fp = fopen("fdm_phase.dat","w");
    for(int iy = 0; iy < ny; iy ++) {
        float lon = lon0 + iy * dy;
        for(int ix = 0; ix < nx; ix ++) {
            float lat = lat0 - ix * dx;
            fprintf(fp,"%f %f %g\n",lon,lat,fdm(ix,iy));
        }
    }
    fclose(fp);

    // group velocity fdm
    for(int iy = 0; iy < ny; iy ++) 
    for(int ix = 0; ix < nx; ix ++) {
        if(iy > ny/3 && iy < ny /3 * 2 && ix > nx / 3 && ix < nx/3*2) {
            vel(ix,iy) = 2.5;
        }
    }
    fmst_reset(vel.data());
    for(int ir = 0; ir < nr; ir ++) {
        fmst_raypath(srcx,srcz,rcx[ir],rcz[ir],fdm.data());
    }
    fmst_finalize();

    fp = fopen("fdm_group.dat","w");
    for(int iy = 0; iy < ny; iy ++) {
        float lon = lon0 + iy * dy;
        for(int ix = 0; ix < nx; ix ++) {
            float lat = lat0 - ix * dx;
            fprintf(fp,"%f %f %g\n",lon,lat,fdm(ix,iy));
        }
    }
    fclose(fp);
}