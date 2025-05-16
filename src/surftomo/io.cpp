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
    std::ifstream infile; infile.open(modfile);
    if(!infile.is_open()) {
        printf("cannot open %s\n",modfile.c_str());
        exit(1);
    }
    std::string line;
    getline(infile,line);
    sscanf(line.c_str(),"%d%d%d",&nx,&ny,&nz);
    getline(infile,line);
    sscanf(line.c_str(),"%f%f",&goxd,&gozd);
    getline(infile,line);
    sscanf(line.c_str(),"%f%f",&dvxd,&dvzd);

    // output some infomation to screen
    printf("\nModel Description:\n");
    printf("===================================\n");
    printf("model origin: latitude,longitude\n");
    printf("   %g   %g\n",goxd,gozd);
    printf("model grid spacing: dlat,dlon\n");
    printf("   %g   %g\n",dvxd,dvzd);
    printf("model dimension: nlat,nlon,nz\n");
    printf("%5d %5d %5d\n",nx,ny,nz); 

    // allocate space
    vsinit.resize(nx,ny,nz);
    fvec depth(nz);
    
    // read depth
    printf("Grid points in depth direction:(km):\n");
    getline(infile,line);
    size_t len = line.size();
    char tmp[len + 10];
    strcpy(tmp,line.c_str());
    char *starp = tmp,*endp = NULL;
    for(int i = 0; i < nz; i ++) {
        depth[i] = std::strtof(starp,&endp);
        starp = endp;
        printf("%7.2f ",depth[i]);
    }
    printf("\n\n");

    // read model
    for(int k=0;k<nz;k++){
    for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++){
        infile >> vsinit(i,j,k);
    }}}
    infile.close();

    // set surftime
    surf.set_model(depth,goxd,gozd,dvxd,dvzd);
    lon.resize(ny); lat.resize(nx); dep = depth;
    for(int i = 0; i < nx; i ++){
        lat[i] = goxd - i * dvxd;
    }
    for(int i = 0; i < ny; i ++){
        lon[i] = gozd + i * dvzd;
    }

    // read true model if required
    if(param.ifsyn && modtrue != "None") {
        vstrue.resize(nx,ny,nz);
        infile.open(modtrue);
        if(!infile.is_open()) {
            printf("cannot open %s\n",modtrue.c_str());
            exit(1);
        }
        getline(infile,line); getline(infile,line); getline(infile,line);
        getline(infile,line);
        for(int k=0;k<nz;k++){
        for(int j=0;j<ny;j++){
        for(int i=0;i<nx;i++){
            infile >> vstrue(i,j,k);
        }}}
        infile.close();
    }
    else {
        printf("You should input a trumodel file (e.g. MOD.true)!\n");
        exit(1);
    }
}

/**
 * @brief read dispersion data
 * 
 * @param datafile 
 */
void DSurfTomo:: 
read_data(const std::string &datafile)
{
    surf.read_swd_data(datafile);
}