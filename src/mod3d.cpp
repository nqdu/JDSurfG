#include"defmod.hpp"
#include"empirical.hpp"
using namespace Eigen;
// mod3d defined
void MOD3d:: gravity(coo_matrix<float> &A,Tensor<float,3> &vsf,VectorXf &dgsyn){
    float rho,a,b; // density, alpha, beta

    if(vs.size()!= vsf.size()){
        std::cout<<"the dimension dismatch!" << std::endl;
        exit(0);
    }
    int size = (ny-2) * (nx-2) * (nz-1);
    VectorXf drho(size);
    dgsyn.setZero();
    for(int k=0;k<nz-1;k++){
        for(int j=0;j< ny-2;j++){
            for(int i=0;i<nx-2;i++){
                b = vsf(i+1,j+1,k);
                int n = k * (ny-2) * (nx-2) + j *(nx-2) + i;
                empirical_relation(&b,&a,&rho);
                drho(n) = rho;
                b = vs(i+1,j+1,k);
                empirical_relation(&b,&a,&rho);
                drho(n) -= rho;
            }
        }
    }
    A.aprod(1,drho.data(),dgsyn.data());
}
