#pragma once
#include"SurfaceWave.hpp"
#include"gravmat.hpp"

class InverseParamsBase
{
    public:
    float smooth,damp; // parameters for lsmr
    float minvel,maxvel; // velocity prior information
    int maxiter; // maxiteration
    int ifsyn; // if or not do checkerboard test
};

class DSurfTomoParams: public InverseParamsBase
{
    public:
    int nthreads;
    float noiselevel; // std of gaussion distribution
};

class JointTomoParams: public InverseParamsBase
{
    public:
    int nthreads;
    float weight1,weight2; // weight parameters for each data
    float noiselevel,noiselevel1; // std of gaussion distribution
    float p; // relative weight params, see Julia(2000)
    bool remove_average;
};

class DSurfTomo{
public:
    int unknowns,num_data;
    MOD3d mod;
    SurfTime surf;
    DSurfTomoParams param;
    Eigen::Tensor<float,3> vstrue;
    
public:
    void readdata(std::string &paramfile,std::string &datafile,
                std::string &modfile,std::string &modture );
    void checkerboard();
    void forward(Eigen::Tensor<float,3> &vs,Eigen::VectorXf &dsyn);
    void inversion(Eigen::Tensor<float,3> &vsf,Eigen::VectorXf &dsyn);
};

class JointTomo
{
public:
    int unknowns,num_data;
    MOD3d mod;
    SurfTime surf;
    Eigen::Tensor<float,3> vstrue;
    OBSSphGraRandom obsg; // gravity observations
    JointTomoParams param; // tomography parameters
    
private:
    MOD3d modref; // reference model
    csr_matrix<float> gmat; // gravity forward computation matrix


    void assemble(Eigen::Tensor<float,3> &vsf,csr_matrix<float> &smat,
                    float weight1,float weight2);

public:
    void readdata(std::string &paramfile,std::string &modfile,std::string &surfdata,
                 std::string &gravdata,std::string &gravmat,std::string &refmod,
                 std::string &modtrue);
    
    void checkerboard();
    void forward(Eigen::Tensor<float,3> &vs,Eigen::VectorXf &dsyn,Eigen::VectorXf &dg);
    void inversion(Eigen::Tensor<float,3> &vsf,Eigen::VectorXf &dsyn,Eigen::VectorXf &dg);
};