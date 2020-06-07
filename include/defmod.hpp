#pragma once
#include<iostream>
#include<string>
#include<unsupported/Eigen/CXX11/Tensor>
#include"coo_matrix.hpp"

// initial model defined
class MOD3d{
    public :
    int nx,ny,nz;
    float goxd,gozd;
    float dvxd,dvzd;
    Eigen::VectorXf lon,lat,dep;
    Eigen::Tensor<float,3> vs;

    void gravity(coo_matrix<float> &A,Eigen::Tensor<float,3> &vsf,Eigen::VectorXf &dgsyn);
};

// surface wave traveltime class
class SurfTime{
    public :
    int kmaxRc,kmaxRg,kmaxLc,kmaxLg; //! number of period points for each wavetype!
    int kmax; // 
    int num_data;
    int n;

    int sublayer;
    Eigen::VectorXd tRc,tRg,tLc,tLg; // period vector
    Eigen::MatrixXf scxf,sczf ;
    Eigen::Tensor<float,3> rcxf,rczf;
    Eigen::MatrixXi wavetype,igrt,nrc1,periods; // wavetype, velotype,no. of receivers
                                         // and periods  per each period and source index 
    Eigen::VectorXi nsrc1; //no. of sources for each period
    Eigen::VectorXf obst,dist; // observed data and distance;

    void forward(MOD3d &mod,Eigen::Tensor<float,3> &vs,Eigen::VectorXf &dsyn);
    void checkerboard(MOD3d &mod,Eigen::Tensor<float,3> &vs,Eigen::VectorXf &dsyn,float noiselevel);
    int FrechetKernel(MOD3d &mod,Eigen::Tensor<float,3> &vsf,
                    Eigen::VectorXf &dsyn,coo_matrix<float> &smat);
    void inversion(MOD3d &mod,Eigen::Tensor<float,3> &vsf,coo_matrix<float> &smat,
                   Eigen::VectorXf &dv,Eigen::VectorXf &dsyn,float damp,float weight,
                   float minvel=-1.0,float maxvel=-1.0);
};

//=== Random observation system, all variables are in Observation-Centred Coordinate System
class OBSSphGraRandom{
    public:
    int   np;     //=== number of observation grids in x, y, directions ===
	float  z0;	//=== observation elevation
    float  *lon, *lat; //=== observation lon/lat

	float *Gr; //=== Observered gravtiy in radial direction.

    int  israd; //=== Flag for lon/lat (ISRAD=0) or lonrad/colatrad/r (ISRAD=1)====

    ~OBSSphGraRandom(){
        delete[] Gr; delete[] lon; delete[] lat;
    }

    void chancoor(int flag);
    void read_obs_data(std::string filename);
};

class TomoParams
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

