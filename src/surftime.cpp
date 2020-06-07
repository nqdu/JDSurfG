#include"defmod.hpp"
#include"calsurf.hpp"
#include"gaussian.hpp"
#include<vector>
using Eigen::Tensor;
using Eigen::VectorXf;
using Eigen::TensorMap;


/*
compute fundamental-mode surface wave traveltimes 
Parameters:
------------------------------------------------
    mod : mod3d class, contains model dimension, and 
          grid information
    vs  : input S-wave velocity model, Tensor class
    dsyn : synthetic traveltimes
*/
void SurfTime :: forward(MOD3d &mod,Tensor<float,3> &vs,VectorXf &dsyn)
{
    int nsrc = scxf.rows();

    synthetic(mod.nx,mod.ny,mod.nz,vs.data(),dsyn.data(),
        mod.goxd,mod.gozd,mod.dvxd,mod.dvzd,kmaxRc,kmaxRg,
        kmaxLc,kmaxLg,tRc.data(),tRg.data(),tLc.data(),tLg.data(),
        wavetype.data(),igrt.data(),periods.data(),mod.dep.data(),
        sublayer,scxf.data(),sczf.data(),rcxf.data(),rczf.data(),
        nrc1.data(),nsrc1.data(),kmax,nsrc,nsrc);
}

/*
synthetic checkerboard traveltimes, synthetic and add noise
Parameters:
------------------------------------------------
    mod  : mod3d class, contains model dimension, and 
          grid information
    vs   : input S-wave velocity model, Tensor class
    dsyn : synthetic traveltimes
    noiselevel : std of gaussian
*/
void SurfTime::checkerboard(MOD3d &mod,Tensor<float,3> &vs,
            VectorXf &dsyn,float noiselevel)
{
    forward(mod,vs,dsyn);
    // add noise
    for(int i=0;i<num_data;i++){
        dsyn(i) += noiselevel * gaussian();
    }

}

/*
compute surface wave frechet kernel and traveltimes
Parameters:
------------------------------------------------
    mod  : mod3d class, contains model dimension, and 
          grid information
    vs   : input S-wave velocity model, Tensor class
    dsyn : synthetic traveltimes
    smat : sparse matrix to store Frechet kernel

*/
int SurfTime ::FrechetKernel(MOD3d &mod,Tensor<float,3> &vs,VectorXf &dsyn,
                             coo_matrix<float> &smat)
{
    int nsrc = scxf.rows();

    // compute frechet kernel
    int nar = 0;
    CalSurfG(mod.nx,mod.ny,mod.nz,smat.nonzeros,vs.data(),&smat.rw[0],&smat.col[0],
            &smat.val[0],dsyn.data(),mod.goxd,mod.gozd,mod.dvxd,mod.dvzd,kmaxRc,
            kmaxRg,kmaxLc,kmaxLg,tRc.data(),tRg.data(),tLc.data(),tLg.data(),
            wavetype.data(),igrt.data(),periods.data(),mod.dep.data(),
            sublayer,scxf.data(),sczf.data(),rcxf.data(),rczf.data(),
            nrc1.data(),nsrc1.data(),kmax,nsrc,nsrc,&nar);
    return nar;
}

/*
Direct Surface Wave Tomography, inversion for one iteration
Parameters:
------------------------------------------------
    mod  : mod3d class, contains model dimension, and 
          grid information
    vs   : input S-wave velocity model, Tensor class
    dsyn : synthetic traveltimes
    smat : sparse matrix to store Frechet kernel and regularization terms
    dv : velocity variation for this iteration

*/
void SurfTime :: inversion(MOD3d &mod,Tensor<float,3> &vs,coo_matrix<float> &smat,
                            VectorXf &dv,VectorXf &dsyn,float damp,float weight,
                            float minvel,float maxvel)
{
    int m = num_data;
    int nx = mod.nx, ny = mod.ny;
    int nz = mod.nz;
    VectorXf res(m + n);

    // compute FrechetKernel
    smat.setZeros();
    int nar = FrechetKernel(mod,vs,dsyn,smat);
    // nar is no. of zeros of current smat, without regularization
    
    // compute residuals
    res.setZero();
    res.segment(0,m) = obst - dsyn;

    // add regularization terms
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
            smat.col[nar] = k * (ny -2) * (nx -2) + j * (nx -2) + i;
            smat.val[nar] = 2.0 * weight;
            smat.rw[nar] = count + m;
            nar ++;
            count ++ ; 
        }
        else{
            if(nar  + 7 > smat.nonzeros){
                std::cout << "please increase sparse ratio!" << std::endl;
                exit(0);
            }
            int rwc = count + m;  // current row
            int clc = k * (ny -2) * (nx -2) + j * (nx -2) + i;// current column
            smat.val[nar] = 6.0 * weight;
            smat.col[nar] = clc;
            smat.rw[nar] = rwc;

            // x direction
            smat.val[nar + 1] = -weight;
            smat.rw[nar + 1] = rwc;
            smat.col[nar + 1] = clc - 1;
            smat.val[nar + 2] = -weight;
            smat.rw[nar + 2] = rwc;
            smat.col[nar + 2] = clc + 1;

            // y direction
            smat.val[nar + 3] = -weight;
            smat.rw[nar + 3] = rwc;
            smat.col[nar + 3] = clc - (nx - 2);
            smat.val[nar + 4] = -weight;
            smat.rw[nar + 4] = rwc;
            smat.col[nar + 4] = clc + (nx - 2);

            // z direction
            smat.val[nar + 5] = -weight;
            smat.rw[nar + 5] = rwc;
            smat.col[nar + 5] = clc - (nx - 2) * (ny - 2);
            smat.val[nar + 6] = -weight;
            smat.rw[nar + 6] = rwc;
            smat.col[nar + 6] = clc + (nx - 2) * (ny - 2);

            nar += 7;
            count++;
        }
    }}}

    // renew sparse matrix meta-informations
    //smat.m = m + count;
    smat.nonzeros = nar;
    smat.cpp2fortran();

    // solve equations by lsmr
    std::cout <<" solving linear systems by LSMR ..." << std::endl;
    int itnlim = n * 2;
    LSMRDict<float> dict(itnlim,10,damp,weight);
    smat.LsmrSolver(res.data(),dv.data(),dict);
    std::cout << "max negative and positive perturbation: " \
                << dv.minCoeff() <<" " << dv.maxCoeff()\
                <<std::endl;
    
    // tackle large variations
    TensorMap<Tensor<float,3>>dx(dv.data(),nx-2,ny-2,nz-1);
    for(int k=0;k<nz-1;k++){
        for(int j=0;j<ny-2;j++){
            for(int i=0;i<nx-2;i++){
                if(dx(i,j,k) > 0.5) dx(i,j,k) = 0.5;
                if(dx(i,j,k) < -0.5) dx(i,j,k) = -0.5;
                float temp = dx(i,j,k) + vs(i+1,j+1,k);
                if(temp > maxvel && maxvel > 0.0) temp = maxvel;
                if(temp < minvel && minvel > 0.0) temp = minvel;
                vs(i+1,j+1,k) = temp;
            }
        }
    }
}
