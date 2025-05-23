#include "surftomo/surftomo.hpp"
#include "shared/IOFunc.hpp"
#include "shared/csr_matrix.hpp"
#include "shared/gaussian.hpp"
#include "shared/smoothing.hpp"

/**
 * @brief inversion for one step by using LSMR method
 * 
 * @param vsf current/next velocity
 * @param dsyn synthetic data for current velocity model
 */
void DSurfTomo::
inversion(fmat3 &vsf,fvec &dsyn)
{
    int m = surf.obst.size();
    int nx = vsf.dimension(0), ny = vsf.dimension(1);
    int nz = vsf.dimension(2);
    int n = (nx-2) * (ny -2) * (nz  - 1);

    // residual vector and velocity variation vector
    fvec res(m + n),dv(n);
    res.setZero();
    dv.setZero();

    // compute and save frechet matrix
    int nar = surf.frechet_matrix(vsf,dsyn,"frechet.bin");
    res.segment(0,m) = surf.obst - dsyn;

    // initialize matrix and read frechet binary file
    csr_matrix smat(m+n,n,nar + n*7);
    std::ifstream fp("frechet.bin",std::ios::binary);
    fp.read((char*)smat.indptr,sizeof(int));
    fp.read((char*)smat.indptr,sizeof(int));
    smat.indptr[0] = 0;
    for(int i=0;i<m;i++){
        int col,nar1;
        fp.read((char*)&col,sizeof(int));
        fp.read((char*)&nar1,sizeof(int));
        int start = smat.indptr[col];
        smat.indptr[col+1] = start +  nar1;

        // read data
        fp.read((char*)(smat.indices + start),sizeof(int) * nar1);
        fp.read((char*)(smat.data + start),sizeof(float) * nar1);
    }
    for(int i=m;i<m+n;i++) smat.indptr[i+1] = smat.indptr[m];

    // delete kernel file
    std::remove("frechet.bin");

    // add regularization
    add_regularization(smat,param.smooth,nx-2,ny-2,nz-1);

    // solve equations by lsmr
    printf("solving linear systems by LSMR ...\n");
    int itnlim = n * 2; 
    LSMRDict<float> dict(itnlim,10,param.damp,param.smooth,param.nthreads);
    dict.verbose = false; 
    smat.lsmr_solver(dv.data(),res.data(),dict);
    printf("max negative and positive perturbation:  %f %f\n",
            dv.minCoeff(),dv.maxCoeff());
    
    // tackle large variations
    Eigen::TensorMap<fmat3>dx(dv.data(),nx-2,ny-2,nz-1);
    float minvel = param.minvel;
    float maxvel = param.maxvel;
    for(int k=0;k<nz-1;k++){
    for(int j=0;j<ny-2;j++){
    for(int i=0;i<nx-2;i++){
        // if(dx(i,j,k) > 0.5) dx(i,j,k) = 0.5;
        // if(dx(i,j,k) < -0.5) dx(i,j,k) = -0.5;
        float temp = dx(i,j,k) + vsf(i+1,j+1,k);
        if(temp > maxvel) temp = maxvel;
        if(temp < minvel) temp = minvel;
        if(vsf(i+1,j+1,k) !=0.0) vsf(i+1,j+1,k) = temp;
    }}}
}

void DSurfTomo :: 
checkerboard()
{   
    int m = surf.obst.size();
    fvec dsyn(m);
    surf.travel_time(vstrue,dsyn);

    // add noise
    for(int i=0;i<m;i++){
        //dsyn(i) *= (1.0 + param.noiselevel * gaussian());
        dsyn(i) += (float) (param.noiselevel * gaussian());
    }
    surf.obst = dsyn*1.0f;
}

float DSurfTomo:: 
compute_misfit(const fvec &dsyn) const
{
    return (dsyn-surf.obst).square().sum() * 0.5;
}

void DSurfTomo:: 
forward(const fvec &x,fvec &dsyn) const
{
    // copy x to vs 
    fmat3 vsf = vsinit;
    int nx = vsf.dimension(0), ny = vsf.dimension(1);
    int nz = vsf.dimension(2);

    int ic = 0;
    for(int k=0;k<nz-1;k++){
    for(int j=0;j<ny-2;j++){
    for(int i=0;i<nx-2;i++){
        vsf(i+1,j+1,k) = x[ic];
        ic += 1;
    }}}

    // forward
    surf.travel_time(vsf,dsyn);
}

void DSurfTomo:: 
compute_grad(const fvec &x,fvec &dsyn,fvec &grad) const
{
    // copy x to vs 
    fmat3 vsf = vsinit;
    int nx = vsf.dimension(0), ny = vsf.dimension(1);
    int nz = vsf.dimension(2);

    int ic = 0;
    for(int k=0;k<nz-1;k++){
    for(int j=0;j<ny-2;j++){
    for(int i=0;i<nx-2;i++){
        vsf(i+1,j+1,k) = x[ic];
        ic += 1;
    }}}

    // compute grad/dsyn
    surf.compute_grad(vsf,dsyn,grad);
}

/**
 * @brief smooth gradient by using gaussian smoothing algorithm
 * 
 * @param grad  gradient array, shape(nlat-2,nlon-2,nz-1), col major
 */
void DSurfTomo ::
smoothing(fvec &grad) const
{
    int nx = vsinit.dimension(0), ny = vsinit.dimension(1);
    int nz = vsinit.dimension(2);

    // find minimum dz
    float dz = 1000 * (dep[nz-2] - dep[0]); 
    for(int i = 0; i < nz-2; i ++) {
        dz = std::min(dz,dep[i+1] - dep[i]);
    }

    // create regular data
    int nz1 = (dep[nz-2] - dep[0]) / (dz * 0.99);
    dz = (dep[nz-2] - dep[0]) / (nz1 - 1);
    fmat3 gradr(nx-2,ny-2,nz1);
    interp_irregular_z(grad.data(),gradr.data(),nx-2,ny-2,nz-1,nz1,dep.data(),true);

    // smoothing
    if (param.smooth_in_km) {
        float dx = std::abs(lat[1] - lat[0]);
        float dy = std::abs(lon[1] - lon[0]);
        smooth_sph_pde(gradr.data(),nx-2,ny-2,nz1,dx,dy,dz,
                        lat[1],lon[1],dep[0],param.sigma_h,
                        param.sigma_v);
    }
    else {
        smooth_cart_pde(gradr.data(),nx-2,ny-2,nz1,param.sigma_h,param.sigma_v);
    }

    // interpolate back
    interp_irregular_z(grad.data(),gradr.data(),nx-2,ny-2,nz-1,nz1,dep.data(),false);
}