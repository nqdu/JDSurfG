#pragma once
#include"defmod.hpp"

// class for direct surface wave tomography
class SurfTomo{
    public:
    int n;
    int num_data;
    SurfTime surf;
    MOD3d mod;
    TomoParams param;

    Eigen::Tensor<float,3> vstrue;

    void readdata(std::string paramfile,std::string datafile,
                std::string modfile,std::string modture );
    

    void checkerboard();
    void forward(Eigen::Tensor<float,3> &vs,Eigen::VectorXf &dsyn);
    void inversion(Eigen::Tensor<float,3> &vsf,Eigen::VectorXf &dv,
                    Eigen::VectorXf &dsyn);
};