#include"bkg_model.hpp"
#include"utils.hpp"
using Eigen::Tensor;
using Eigen::VectorXf;

/** Compute gravity anomalies for a given velocity model
 * Parameters:
 * ---------------------------------------------
 *  A       : csr_matrix<float>, gravity matrix
 * vsf      : velocity model
 * dgsyn    : synthetic gravity data
 */
void MOD3d:: gravity(csr_matrix<float> &A,Tensor<float,3> &vsf,VectorXf &dgsyn){
    float rho,a,b; // density, alpha, beta

    if(vs.size()!= vsf.size()){
        std::cout<<"the dimension dismatch!" << std::endl;
        exit(0);
    }
    int size = (ny-2) * (nx-2) * (nz-1);
    VectorXf drho(size);
    dgsyn.setZero();
    for(int k=0;k<nz-1;k++){
    for(int j=0;j<ny-2;j++){
    for(int i=0;i<nx-2;i++){
        b = vsf(i+1,j+1,k);
        int n = k * (ny-2) * (nx-2) + j *(nx-2) + i;
        empirical_relation(&b,&a,&rho);
        drho(n) = rho;
        b = vs(i+1,j+1,k);
        empirical_relation(&b,&a,&rho);
        drho(n) -= rho;
    }}}
    A.aprod(1,drho.data(),dgsyn.data());
}

/**
 * Add 2-th order Tikhonov Regularization Matrix to a csr_matrix
 * Parameters:
 * ----------------------------------------------
 * smat     : csr_matrix<float>
 * weight   : smooth factor
 * 
 * Examples:
 * ----------------------------------------------------
 * For 1-D problem with dimension n = 5, regularization matrix is:
 *             [ 2   0  0  0  0 ]
 *             [ -1  2 -1  0  0 ]
 *             [ 0   -1 2  -1 0 ]
 *             [ 0   0  -1  2 -1]
 *             [ 0   0   0  0  2]
 */
void MOD3d::add_regularization(csr_matrix<float> &smat,float weight)
{
    // get parameters required
    int n = (nx-2) * (ny -2) * (nz  - 1); // model dimension
    int nar = smat.nonzeros - n * 7; // nonzeros excluding smooth term
    int m = smat.rows() - n; // data dimension

    int count = 0;
    for(int k=0;k<nz-1;k++){
    for(int j=0;j<ny-2;j++){
    for(int i=0;i<nx-2;i++){
        if( i==0 || i==nx-3 || j==0 || j==ny-3 || k==0 || k==nz-2){
            // and more restrictions to boundary points
            if(nar + 1 > smat.nonzeros){
                std::cout << "please increase sparse ratio!" << std::endl;
                exit(0);
            }
            int rwc = count + m;
            smat.data[nar] = 2.0 * weight;
            smat.indices[nar] = k * (ny -2) * (nx -2) + j * (nx -2) + i;
            smat.indptr[rwc + 1] = smat.indptr[rwc] + 1;
            nar += 1;
            count += 1;
        }
        else{
            if(nar  + 7 > smat.nonzeros){
                std::cout << "please increase sparse ratio!" << std::endl;
                exit(0);
            }
            int rwc = count + m;  // current row
            smat.indptr[rwc +1] = smat.indptr[rwc] + 7;
            int clc = k * (ny -2) * (nx -2) + j * (nx -2) + i;// current column
            smat.data[nar] = 6.0 * weight;
            smat.indices[nar] = clc;

            // x direction
            smat.data[nar + 1] = -weight;
            smat.indices[nar + 1] = clc - 1;
            smat.data[nar + 2] = -weight;
            smat.indices[nar + 2] = clc + 1;

            // y direction
            smat.data[nar + 3] = -weight;
            smat.indices[nar + 3] = clc - (nx - 2);
            smat.data[nar + 4] = -weight;
            smat.indices[nar + 4] = clc + (nx - 2);

            // z direction
            smat.data[nar + 5] = -weight;
            smat.indices[nar + 5] = clc - (nx - 2) * (ny - 2);
            smat.data[nar + 6] = -weight;
            smat.indices[nar + 6] = clc + (nx - 2) * (ny - 2);

            nar += 7;
            count += 1;
        }
    }}}
    smat.nonzeros = nar;
}

void MOD3d::
empirical(float plat,float plon,float pdep,float beta,
                            float &alpha,float &rho)
{

}