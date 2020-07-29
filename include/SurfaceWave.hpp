#pragma once
#include<iostream>
#include<unsupported/Eigen/CXX11/Tensor>
#include<string>
#include"bkg_model.hpp"

class StationPair
{
    public:
    std::string wavetype;
    int counter; // data counter
    int period_idx; // which period, 0 for minimum period used and >0 for else
    //int mode_indx; // which mode, 0 for fundamental and >=1 for higher modes
    int nr; // no. of receivers

    public:
    float srcx,srcz; // source station colat and lon, in rad
    std::vector<float> rcx,rcz;// receiver station coordinates, colat and lon ,in rad
    std::vector<float> dist,obstime; // distance and traveltime

    StationPair(std::string wavetp,int count,int p_idx,int nrcv,float scx,float scz)
    {
        wavetype = wavetp;
        counter = count;
        period_idx = p_idx;
        nr = nrcv;
        srcx = scx;
        srcz = scz;
        rcx.resize(nr);
        rcz.resize(nr);
        dist.resize(nr);
        obstime.resize(nr);
    }

};

class SurfTime
{
    public:
    int num_data; 
    int kmaxRc,kmaxRg,kmaxLc,kmaxLg; // no. of periods used in each mode
    int kmax;
    Eigen::VectorXd tRc,tRg,tLc,tLg; // period vector for each mode
    Eigen::VectorXf obst; // observations
    Eigen::VectorXf sta_dist;

    public:
    int sublayer;
    std::vector<StationPair> Pairs;

    private:
    void DipersionMap(MOD3d &mod,Eigen::Tensor<float,3> &vs,Eigen::MatrixXd &pv);
    void SurfWaveKernel(MOD3d &mod,Eigen::Tensor<float,3> &vs,Eigen::MatrixXd &pv,
                           Eigen::Tensor<double,3> &kernel );
    int get_period_index(int idx,std::string wavetype);

    public:
    int read_data(std::string filename);
    void TravelTime(MOD3d &mod,Eigen::Tensor<float,3> &vs,Eigen::VectorXf &data);
    int FrechetKernel(MOD3d &mod,Eigen::Tensor<float,3> &vs,Eigen::VectorXf &data,
                        std::string save_path);
};
