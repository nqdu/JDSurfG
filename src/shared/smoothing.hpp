#ifndef JDSURFG_SHARED_SMOOTHING_H_
#define JDSURFG_SHARED_SMOOTHING_H_

void smooth_pde_irregular(
    float* __restrict gradorg,
    int nx,int ny,int nz, const float *depth_ireg,
    float lat0,float lon0,
    float dx,float dy,
    float sigma_h,float sigma_v,
    bool smooth_in_km
);

#endif // end JDSURFG_SHARED_SMOOTHING_H_
