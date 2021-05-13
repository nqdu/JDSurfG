#pragma once
#include"SurfaceWave.hpp"
#include"gravmat.hpp"
class DSurfTomoParams
{
    public:
    float smooth,damp; // parameters for lsmr
    float spra;  // sparse matrix sparse ratio
    int maxiter; // maxiteration
    float minvel,maxvel; // velocity prior information
    int ifsyn; // if or not do checkerboard test
    float noiselevel; // std of gaussion distribution
};

class JointTomoParams
{
    public:
    float smooth,damp; // parameters for lsmr
    float weight1,weight2; // weight parameters for each data
    float spra;  // sparse matrix sparse ratio
    int maxiter; // maxiteration
    float minvel,maxvel; // velocity prior information
    int ifsyn; // if or not do checkerboard test
    float noiselevel,noiselevel1; // std of gaussion distribution
    float p; // relative weight params, see Julia(2000)
};

class DSurfTomo{
    public:
    MOD3d mod;
    SurfTime surf;
    DSurfTomoParams param;
    Eigen::Tensor<float,3> vstrue;
    int unknowns,num_data;
    
    public:
    int readdata(std::string &paramfile,std::string &datafile,
                std::string &modfile,std::string &modture );
    void checkerboard();
    void forward(Eigen::Tensor<float,3> &vs,Eigen::VectorXf &dsyn);
    void inversion(Eigen::Tensor<float,3> &vsf,Eigen::VectorXf &dsyn);
};

class JointTomo{
    public:
    int num_data;
    int n;

    SurfTime surf; // surface wave traveltime class
    OBSSphGraRandom obsg; // gravity observations
    MOD3d mod; // initial model class
    Eigen::Tensor<float,3> vstrue;
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