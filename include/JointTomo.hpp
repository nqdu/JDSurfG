#pragma once
#include"defmod.hpp"

class JointTomo{
    public:
    int num_data;
    int n;

    SurfTime surf; // surface wave traveltime class
    MOD3d mod; // initial model class
    MOD3d modref; // reference model
    TomoParams param; // tomography parameters
    OBSSphGraRandom obsg; // gravity observations
    coo_matrix<float> gmat; // gravity forward computation matrix

    Eigen::Tensor<float,3> vstrue;

    void readdata(std::string paramfile,std::string modfile,std::string surfdata,
                 std::string gravdata,std::string gravmat,std::string refmod,
                 std::string modtrue);
    
    void checkerboard();
    void forward(Eigen::Tensor<float,3> &vs,Eigen::VectorXf &dsyn,Eigen::VectorXf &dg);
    void inversion(Eigen::Tensor<float,3> &vsf,Eigen::VectorXf &dv,
                    Eigen::VectorXf &dsyn,Eigen::VectorXf &dg);

};