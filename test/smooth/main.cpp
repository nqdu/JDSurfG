#include "numerical.hpp"
#include "shared/smoothing.hpp"
#include "shared/spherical.hpp"
#include <iostream>
#include <chrono>

int main () {
    int nx,ny,nz;
    float dx,dy,dz,lat0,lon0,z0;
    float sigma_h,sigma_v;

    // read parameters
    FILE *fp = fopen("smooth.in","r");
    char line[256];
    fgets(line,sizeof(line),fp);  sscanf(line,"%d%d%d",&nx,&ny,&nz);
    fgets(line,sizeof(line),fp);  sscanf(line,"%f%f%f",&dx,&dy,&dz);
    fgets(line,sizeof(line),fp);  sscanf(line,"%f%f%f",&lat0,&lon0,&z0);
    fgets(line,sizeof(line),fp);  sscanf(line,"%f%f",&sigma_h,&sigma_v);
    fclose(fp);

    // gradient 
    fmat3 grad(nx,ny,nz); grad.setConstant(-1.);
    for(int iz = 0; iz < nz; iz ++)
    for(int iy = 0; iy < ny; iy ++)
    for(int ix = 0; ix < nx; ix ++) {
       if(ix > nx /3 && ix < nx/3 * 2 && iy > ny / 3 && iy < ny / 3 * 2 && iz > nz/3 && iz < nz/3*2) {
            grad(ix,iy,iz) = 1.;
        }
    }

    // save function
    auto save_grad = [nz,ny,nx,lat0,lon0,z0,dx,dy,dz]
                (const char *filename,fmat3 &grad,int nx0,int nx1,
                int ny0,int ny1,int nz0,int nz1) 
    {
        const float deg2rad = M_PI / 180.;
        FILE *fp = fopen(filename,"w");
        for(int iz = nz0; iz < nz1; iz ++) {
            float z = z0 + iz * dz;
            for(int iy = ny0; iy < ny1; iy ++) {
                float lon = lon0 + dy * iy;
                for(int ix = nx0; ix < nx1; ix ++) {
                    float lat = lat0 - dx * ix;

                    float d = delsph((90-lat) * deg2rad,deg2rad*lon,(90-lat0) * deg2rad,lat0*deg2rad);
                    //d = lat;

                    fprintf(fp,"%f %f %f\n",d,z,grad(ix,iy,iz));
                }
            }
        }
        fclose(fp);
    };

    // save gradient slice 
    save_grad("grad.org.dat",grad,0,nx,ny/2,ny/2+1,0,nz);

    // smooth by using different method
    fmat3 grad1 = grad;
    auto start = std::chrono::steady_clock::now();
    smooth_sph(grad1.data(),nx,ny,nz,dx,dy,dz,lat0,lon0,z0,sigma_h,sigma_v);
    auto end = std::chrono::steady_clock::now();
    float t = std::chrono::duration<float>(end-start).count();
    printf("convolve smoothing time elapsed = %f \n",t);

    start = std::chrono::steady_clock::now();
    smooth_sph_pde(grad.data(),nx,ny,nz,dx,dy,dz,lat0,lon0,z0,sigma_h,sigma_v);
    end = std::chrono::steady_clock::now();
    t = std::chrono::duration<float>(end-start).count();
    printf("PDE smoothing time elapsed = %f \n",t);

    // save file 
    save_grad("grad.smooth.pde.dat",grad,0,nx,ny/2,ny/2+1,0,nz);
    save_grad("grad.smooth.conv.dat",grad1,0,nx,ny/2,ny/2+1,0,nz);
}