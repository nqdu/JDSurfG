#ifndef JDSURFG_SWD_SURFACEWAVE_H_
#define JDSURFG_SWD_SURFACEWAVE_H_

#include "numerical.hpp"
#include <string>

// FMST extensions
extern "C"{
/**
 * @brief initialize fmst lib 
 * @param nx,ny grid dimension for lat/lon direction
 * @param goxdf,gozdf upper left points lat/lon, in deg 
 * @param dxvdf,dvzdf lat/lon interval, in deg 
*/
void fmst_init(int nx,int ny,float goxdf,float gozdf,
            float dvxdf,float dvzdf);

/**
 * @brief destroy memory in fmst lib
*/
void fmst_finalize();

/**
 * @brief compute travel time field 
 * @param velf input velocity model, shape(nlat,nlon), column major
 * @param srcx/z source colat/lon, in rad
 * @param rcx/z receiver source colat/lon, in rad, shape(nr)
 * @param nr # or receivers
 * @param ttime # travel time at receivers, shape(nr)
*/
void fmst_run(const double *velf,float srcx,float srcz,
                const float *rcx,const float *rcz,int nr,float *ttime);

/**
 * @brief reset velocity model as velin
 * 
 * @param velin 
 */
void fmst_reset(const double *velin);

/**
 * @brief compute traveltime/frechet kernel from a given s-r setting by using raypath
 * @param velf input velocity model, shape(nlat,nlon), column major, maybe differ from fmst_travel
 * @param srcx/z source colat/lon, in rad
 * @param rcx/z receiver source colat/lon, in rad, shape(nr)
 * @param fdm # frechet kernel, shape(nlat,nlon), col major
 * @return t travel time for this s-r setting
*/
float fmst_raypath(float scx,float scz,float rcx,float rcz,float *fdm);

}

struct StationPair
{
public:
    std::string wavetype;
    int counter; // data counter
    int period_idx; // which period, 0 for minimum period used and >0 for else
    int mode; // which mode, 0 for fundamental and >=1 for higher modes
    int nr; // no. of receivers

public:
    float srcx,srcz; // source station colat and lon, in rad
    std::vector<float> rcx,rcz;// receiver station coordinates, colat and lon ,in rad
    std::vector<float> dist,obstime; // distance and traveltime

    StationPair(std::string wavetp,int count,int mode,int p_idx,int nrcv,float scx,float scz)
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
        this -> mode = mode;
    }

};

class SurfTime
{
public:
    ivec kmaxRc,kmaxRg,kmaxLc,kmaxLg;
    int kmax;
    dmat2 tRc,tRg,tLc,tLg; // period vector for each mode tRc(mode,nt)
    fvec obst,sta_dist; // observations and gc distance in km
    std::vector<StationPair> Pairs;

private:
    // depth vector for compute dispersion
    fvec depth;

    // model parameters
    float lon0,lat0; // upper left lon/lat in deg
    float dlon,dlat; // model spacing, in deg

private:
    void get_2d_map(const fmat3 &vs,dmat2 &vc,dmat2 &vout) const;
    void get_1d_kernel(const fmat3 &vs,dmat2 &vc,dmat2 &vout,dmat3 &kernel ) const;
    int get_period_index(int idx,int mode,const std::string &wavetp) const ;
    // void SurfWaveKernel(MOD3d &mod,const fmat3 &vs,Eigen::MatrixXd &pv,
    //                        Eigen::Tensor<double,3> &kernel );
    // int get_period_index(int idx,std::string &wavetype);

public:
    void travel_time(const fmat3 &vs,fvec &data) const;
    int frechet_matrix(const fmat3 &vs,fvec &data,const std::string &outfile) const;
    void compute_grad(const fmat3 &vs,fvec &data,fvec &grad) const;
    void write_syn(const fvec &dsyn, const std::string &outfile) const;
    void set_model(const fvec &dep,float goxd,float gozd,float dvxd,float dvzd);
    void read_swd_data(const std::string &datafile);
};

#endif // end JDSURFG_SWD_SURFACEWAVE_H_
