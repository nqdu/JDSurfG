#include <iostream>

#include "surftomo/surftomo.hpp"
#include "shared/IOFunc.hpp"


/**
 * @brief read inverse parameters
 * 
 * @param paramfile parameter file
 */
void DSurfTomo::
read_invparams(const std::string &paramfile)
{
    param.read_file(paramfile);

    // read noise level 
    std::ifstream infile; infile.open(paramfile);
    read_par_regex("NOISE_LEVEL",param.noiselevel,infile);
    infile.close();
}

/**
 * @brief read model from model file
 * 
 * @param modfile initial model file
 * @param modtrue true model file, if = NONE, use the average of initial one
 */
void DSurfTomo ::
read_model(const std::string &modfile,const std::string &modtrue)
{
    int nx,ny,nz;
    float goxd,gozd,dvxd,dvzd;
    fvec depth;
    read_velocity_model(modfile,vsinit,depth,lon,lat,true);
    nx = vsinit.dimension(0);
    ny = vsinit.dimension(1);
    nz = vsinit.dimension(2);

    // read true model if required
    if(param.ifsyn) {
        if(modtrue != "None") {
            fvec depth1,lat1,lon1;
            read_velocity_model(modtrue,vstrue,depth1,lon1,lat1,false);
            
            // check if vstrue has the same size as vsinit
            if(vstrue.dimensions() != vsinit.dimensions() || 
               depth1.size() != depth.size() ||
               lat1.size() != lat.size() ||
               lon1.size() != lon.size()) {
                printf("true model coverage not match initial model coverage!\n");
                printf("please check the true model file: %s\n",modtrue.c_str()); 
                printf("lon size: %td %td\n",lon1.size(),lon.size());
                printf("lat size: %td %td\n",lat1.size(),lat.size());
                printf("depth size: %td %td\n",depth1.size(),depth.size());
                exit(1);
            }
        }
        else {
            printf("You should input a trumodel file (e.g. MOD.true)! when enabling SYN_TEST \n");
            exit(1);
        }
    }

    // set depth 
    dep.resize(nx,ny,nz);        
    for(int k=0;k<nz;k++){
    for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++){
        dep(i,j,k) = depth[k];
    }}}

    // read topography if required
    if(param.topo_corr == 1) {
        printf("\ntopography correction is applied.\n");
        printf("reading topography from topography.dat\n");

        std::string topofile = "topography.dat";
        interpolate_topo(topofile,lat,lon,dep);
    }

    // set 1-D model
    swsol.set_model(dep,lat[0],lon[0],lat[0]-lat[1],lon[1]-lon[0]);
}

/**
 * @brief read dispersion data
 * 
 * @param datafile 
 */
void DSurfTomo:: 
read_data(const std::string &datafile)
{
    swsol.read_swd_data(datafile);
}