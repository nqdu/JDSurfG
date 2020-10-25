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

    /*
    MOD3d & operator = (const MOD3d &md){
        nx = md.nx;
        ny = md.ny;
        nz = md.nz; 
        goxd = md.goxd; gozd = md.gozd;
        dvxd = md.dvxd; dvzd = md.dvzd;
        lon = md.lon; lat = md.lat;
        dep = md.dep;
        vs = md.vs;

        return *this;
    }*/

    void add_regularization(csr_matrix<float> &smat,float smooth);
};