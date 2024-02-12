#include "numerical.hpp"
#include "shared/spherical.hpp"
#include <Eigen/LU>

void smooth_grad(float* grad,int nx,int ny,int nz,float sigma_h,float sigma_v)
{
    float grad1[nx*ny*nz];
    memcpy(grad1,grad,sizeof(float)*nx*ny*nz); 

    //omp_set_num_threads(3);
    #pragma omp parallel for shared(grad,grad1) collapse(3)
    for(int k=0;k<nz;k++){
    for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++){
        int n = k * ny * nx + j * nx + i;
        float sum_k = 0., sum_n = 0.;
        for(int k1=0;k1<nz;k1++){
        for(int j1=0;j1<ny;j1++){
        for(int i1=0;i1<nx;i1++){
            int n1 = k1 * ny * nx + j1 * nx + i1;
            // compute azimuthal distance
            float delta = std::hypot(1. * (i-i1),1. * (j-j1));

            // compute integrals
            float tmp =  std::exp(- pow(delta,2) * 0.5 / (sigma_h * sigma_h))
                         * std::exp(-pow(k-k1,2) * 0.5 / (sigma_v * sigma_v));
            sum_k += grad1[n1] *tmp;
            sum_n += tmp;
        }}}
        grad[n] = sum_k / sum_n;
    }}}
}

void smooth_sph(float* grad,int nx,int ny,int nz,
                float dlat,float dlon,float dz,
                float lat0,float lon0,float z0,
                float sigma_h,float sigma_v)
{
    float grad1[nx*ny*nz];
    memcpy(grad1,grad,sizeof(float)*nx*ny*nz); 
    const float R = 6371.;
    const float deg2rad = M_PI / 180.;

    //omp_set_num_threads(3);
    #pragma omp parallel for shared(grad,grad1) collapse(3)
    for(int k=0;k<nz;k++){
    for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++){
        //coordinates
        float colat = deg2rad * (90. - lat0 + i * dlat);
        float lon = deg2rad * (lon0 + j * dlon);
        float r =  R - (z0 + k * dz);

        // index
        int n = k * ny * nx + j * nx + i;
        float sum_k = 0., sum_n = 0.;
        for(int k1=0;k1<nz;k1++){
        for(int j1=0;j1<ny;j1++){
        for(int i1=0;i1<nx;i1++){
            int n1 = k1 * ny * nx + j1 * nx + i1;

            float colat1 = deg2rad * (90. - lat0 + i1 * dlat);
            float lon1 = deg2rad * (lon0 + j1 * dlon);
            float r1 =  R - (z0 + k1 * dz);

            // compute azimuthal distance
            float delta = delsph(colat,lon,colat1,lon1) / R ;

            // compute integrals
            float tmp =  std::exp(- pow(delta * r1,2) * 0.5 / (sigma_h * sigma_h))
                         * std::exp(-pow(r-r1,2) * 0.5 / (sigma_v * sigma_v));
            sum_k += grad1[n1] *tmp;
            sum_n += tmp;
        }}}
        grad[n] = sum_k / sum_n;
    }}}
}


static int factorial(int n)
{
    int s = 1;
    if(n==0 || n ==1) return s;
    for(int i=2;i<n+1;i++){
        s = s * i;
    } 

    return s ;
}

/**
 * compute central FD coefficients for a specific accuracy order 
 * @param acc order of accuarcy 
 * @param deriv derivative order 
 * @return coef 
*/
void fdcoefs(int acc,float* __restrict__ coefs,int deriv)
{
    // make sure a even number
    int norder = acc; 
    int d = deriv;
    int n = norder / 2;
    assert(n * 2 == norder);

    // allocate space 
    Eigen::MatrixXf A(norder+1,norder+1);
    Eigen::VectorXf b(norder + 1);
    for(int i=0;i<norder+1;i++){
    for(int j=0;j<norder+1;j++){
        A(i,j) = (float) std::pow(j-n,i) / factorial(i);
    }}
    b(d) = 1.0;

    // solve linear system
    Eigen::VectorXf x = A.partialPivLu().solve(b);

    if (d/2 * 2 != d) x(n) = 0.0;

    memcpy(coefs,x.data(),sizeof(float)*(acc+1));
}

/**
 * @brief gaussian smoothing of a 3D array, in cartesian coordinates, the grid interval is 1
 * @param gradin input array, shape(nx,ny,nz), col major
 * @param nx,ny,nz x/y/z nodes number
 * @param sigma_h/v smooth parameter in x/y and z direction
*/
void smooth_cart_pde(float * __restrict__ gradin,int nx,int ny,int nz,float sigma_h,float sigma_v)
{
    Eigen::TensorMap<fmat3> grad(gradin,nx,ny,nz);
    
    // FD coefs
    int norder = 8 ;
    float *coefs = new float[norder+1];
    fdcoefs(norder,coefs,2);
    int n2 = norder / 2;
    int nx1 = nx + n2 * 2, ny1 = ny + n2 * 2, nz1 = nz + n2 * 2;

    // lambda for index conversion
    auto func = [](int i,int nx) {
        int i1 = i;
        if(i1 < 0) i1 = 0;
        if(i1 >=nx) i1 = nx-1;

        return i1;
    };

    // allocate space for a larger mesh
    fmat3 grad_old(nx1,ny1,nz1); grad_old.setZero();
    for(int k=-n2;k<nz+n2;k++){
    for(int j=-n2;j<ny+n2;j++){
    for(int i=-n2;i<nx+n2;i++){
        int k1 = func(k,nz);
        int j1 = func(j,ny);
        int i1 = func(i,nx);
        grad_old(i+n2,j+n2,k+n2) = grad(i1,j1,k1);
    }}}
    fmat3 grad_new = grad_old;

    // compute FD nteps
    const float eta = 1. / 12;
    float sh2 = sigma_h * sigma_h, sv2 = sigma_v * sigma_v;
    int nstep = std::max(sh2,sv2) / (2. * eta);

    // compute new version of sigma_h/v
    float Ch = sh2 / std::max(sh2,sv2),Cv = sv2 / std::max(sh2,sv2);
    Ch *= eta; Cv *= eta;

    for(int it = 0; it < nstep; it ++) {
        for(int k=0;k<nz1;k++){
        for(int j=0;j<ny1;j++){
        for(int i=0;i<nx1;i++){
            bool flag = (i >= n2) && (i < nx1 - n2) &&
                        (j >= n2) && (j < ny1 - n2) &&
                        (k >= n2) && (k < nz1 - n2);
            if(!flag) continue;
            double tmpx = 0., tmpy = 0., tmpz = 0.;
            for(int l = -n2; l <= n2; l ++) {
                tmpx += coefs[l+n2] * grad_old(i+l,j,k);
                tmpy += coefs[l+n2] * grad_old(i,j+l,k);
                tmpz += coefs[l+n2] * grad_old(i,j,k+l);
            }
            grad_new(i,j,k) = grad_old(i,j,k) + Ch * tmpx + Ch * tmpy + Cv * tmpz;
        }}}
        grad_old = grad_new;
    }

    // copy back
    for(int k=0;k<nz;k++){
    for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++){
        grad(i,j,k) = grad_old(i+n2,j+n2,k+n2);
    }}}

    // delete
    delete[] coefs;
}

/**
 * @brief gaussian smoothing of a 3D array, in spherical coordinates
 * @param gradin input array, shape(nx,ny,nz), col major
 * @param nx,ny,nz lat/lon/r nodes number
 * @param lat0/lon0/z0 upper left lat/lon and minimum depth, in deg/km
 * @param dx/dy/dz grid interval in each direction, in deg/km
 * @param sigma_h/v smooth parameter in theta/phi and r direction
*/
void smooth_sph_pde(float* gradin,int nx,int ny,int nz,
                float dx,float dy,float dz,
                float lat0,float lon0,float z0,
                float sigma_h,float sigma_v)
{
    const float R = 6371.;
    const float deg2rad = M_PI / 180.;
    Eigen::TensorMap<fmat3> grad(gradin,nx,ny,nz);

    // FD coefs
    const int norder = 6 ;
    float coefs1[norder + 1],coefs2[norder+1];
    fdcoefs(norder,coefs2,2);
    fdcoefs(norder,coefs1,1);
    const int n2 = norder / 2;
    int nx1 = nx + n2 * 2, ny1 = ny + n2 * 2, nz1 = nz + n2 * 2;

    // lambda for index conversion
    auto func = [](int i,int nx) {
        int i1 = i;
        if(i1 < 0) i1 = 0;
        if(i1 >=nx) i1 = nx-1;

        return i1;
    };

    // allocate space for a larger mesh
    fmat3 grad_old(nx1,ny1,nz1); grad_old.setZero();
    for(int k=-n2;k<nz+n2;k++){
    for(int j=-n2;j<ny+n2;j++){
    for(int i=-n2;i<nx+n2;i++){
        int k1 = func(k,nz);
        int j1 = func(j,ny);
        int i1 = func(i,nx);
        grad_old(i+n2,j+n2,k+n2) = grad(i1,j1,k1);
    }}}
    fmat3 grad_new = grad_old;

    // compute useful arrays
    fvec radi(nz1), theta(nx1),phi(ny1);
    for(int iz = -n2; iz < nz+n2; iz ++) radi[iz+n2] = R - (iz * dz + z0);
    for(int iy = -n2; iy < ny+n2; iy ++) phi[iy+n2] = deg2rad * (lon0 + (iy * dy));
    for(int ix = -n2; ix < nx+n2; ix ++) theta[ix+n2] = deg2rad * (90. - lat0 + (ix * dx));
    float dr = radi[1] - radi[0], dth = theta[1] - theta[0], dphi = phi[1] - phi[0];

    // compute minimal interval in FD mesh
    float h2 = R * 1000;
    for(int k=0;k<nz1;k++){
    for(int j=0;j<ny1;j++){
    for(int i=0;i<nx1;i++){
        float d1 = std::abs(dr);
        float d2 = radi[k] * (dth);
        float d3 = radi[k] * sin(theta[i]) * dphi;
        h2 = std::min(h2,d1);
        h2 = std::min(h2,d2);
        h2 = std::min(h2,d3);
    }}}
    h2 = h2 * h2;

    // compute FD nteps
    const float eta = 1. / 7;
    float sh2 = sigma_h * sigma_h, sv2 = sigma_v * sigma_v;
    int nstep = std::max(sh2,sv2) / (2. * eta * h2);

    // compute new version of sigma_h/v
    float Ch = sh2 / std::max(sh2,sv2),Cv = sv2 / std::max(sh2,sv2);
    Ch *= eta * h2; Cv *= eta * h2;

    for(int it = 0; it < nstep; it ++) {
        // if(it == 0 || it < 10 || it == nstep-1 || (it+1) % 100 == 0){
        //     Eigen::Map<fvec> v(grad_new.data(),grad_new.size());
        //     float per = (it+1.) / nstep * 100;
        //     printf("%d of %d: %g%%, vmin = %g vmax =%g\n",it+1,nstep,per,
        //           v.minCoeff(),v.maxCoeff());
        // }
        #pragma omp parallel for collapse(3)
        for(int k=0;k<nz1;k++){
        for(int j=0;j<ny1;j++){
        for(int i=0;i<nx1;i++){
            bool flag = (i >= n2) && (i < nx1 - n2) &&
                        (j >= n2) && (j < ny1 - n2) &&
                        (k >= n2) && (k < nz1 - n2);
            if(!flag) continue;
            double tmpx1 = 0., tmpx2 = 0., tmpy = 0., tmpz = 0.;
            for(int l = -n2; l <= n2; l ++) {
                tmpx2 += coefs2[l+n2] * grad_old(i+l,j,k); // theta d^2f / dth^2
                tmpx1 += coefs1[l+n2] * grad_old(i+l,j,k); // theta df / dth
                tmpy += coefs2[l+n2] * grad_old(i,j+l,k); // phi df/dphi
                tmpz += coefs2[l+n2] * radi[k+l] * grad_old(i,j,k+l); // r
            }

            grad_new(i,j,k) = grad_old(i,j,k) + 
                              Ch * tmpx2 / (radi[k] * radi[k]) / (dth * dth)+ 
                              Ch * tmpx1 / dth * cos(theta[i]) / (radi[k] * radi[k] * sin(theta[i])) + 
                              Ch * tmpy / (dphi * dphi) / std::pow(radi[k] * sin(theta[i]),2) +
                              Cv * tmpz /radi[k] / (dr*dr) ;
        }}}
        grad_old = grad_new;
    }

    // copy back
    for(int k=0;k<nz;k++){
    for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++){
        grad(i,j,k) = grad_old(i+n2,j+n2,k+n2);
    }}}
}

// void find_quadratic_min(const float *x,const float *y,int n,float &xmin,float &ymin)
// {
//     Eigen::MatrixXd A(3,3);
//     Eigen::VectorXd b(3);
//     A.setConstant(1.);
//     for(int i =0; i < 3; i ++) {
//         A(i,0) = x[i] * x[i];
//         A(i,1) = x[i];
//         b[i]  = y[i];
//     }

//     // solver linear system
//     Eigen::Vector3d m = A.colPivHouseholderQr().solve(b);

//     xmin =  -m(1) / (2.0 * m(0)); // x0 = -b / (2 * a)
//     ymin = (4 * m(0) * m(2) - m(1)*m(1)) / (4 * m(0));

//     if(xmin < 0){
//         int ic = 0;
//         for(int i = 1; i < n; i ++) {
//             if(y[ic] > y[i]) ic = i;
//         }
//         xmin = x[ic];
//         ymin = y[ic];
//     }
// }
