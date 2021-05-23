#pragma once 
#include<unsupported/Eigen/CXX11/Tensor>
#include"csr_matrix.hpp"

// initial model defined
class MOD3d{
    public :
    int nx,ny,nz;
    float goxd,gozd;
    float dvxd,dvzd;
    Eigen::VectorXf lon,lat,dep;
    Eigen::Tensor<float,3> vs;

    void gravity(csr_matrix<float> &A,Eigen::Tensor<float,3> &vsf,Eigen::VectorXf &dgsyn);
    void add_regularization(csr_matrix<float> &smat,float smooth);

    // empirical relations
    void empirical_relation(float vsz,float &vpz,float &rhoz);
    void empirical_deriv(float vp,float vs,float &drda,float &dadb);
};