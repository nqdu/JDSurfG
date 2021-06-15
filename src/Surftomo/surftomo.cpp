#include"tomography.hpp"
#include"IOFunction.hpp"
#include"utils.hpp"
#include"csr_matrix.hpp"
#include<fstream>
using Eigen::Tensor;
using Eigen::VectorXf;
using Eigen::VectorXi;

/**
 * read all input parameters, observations and initial model
 * @param paramfile file contains parameter
 * @param datafile file contains observed dispersion data
 * @param modfile initial model
 * @param modtrue True model, for synthetic test 
 */
void DSurfTomo:: readdata(std::string &paramfile,std::string &datafile,
                        std::string &modfile,std::string &modtrue )
{
    // read paramfile
    std::ifstream infile;
    std::string line;
    std ::istringstream info; 
    infile.open(paramfile);

    // read all parameters
    param.nthreads =  __read_parfile(infile,paramfile,mod,surf,param);
    skipread(infile,line,"%f",&param.noiselevel);
    infile.close();

    // allocate data parameters
    int nx = mod.nx,ny = mod.ny, nz= mod.nz;
    unknowns = (nx-2) * (ny -2) * (nz -1 );
    if(param.ifsyn) vstrue.resize(nx,ny,nz);
    
    // read surface wave traveltime data
    __read_traveltime(datafile,surf);
    this->num_data = surf.num_data;

    // read initial model
    __read_InputModel(modfile,mod.vs,mod.dep.data(),true);

    // read ture model if required
    if(param.ifsyn == 1){
        float z[nz];
        __read_InputModel(modtrue,vstrue,z);
    }

    std::cout << std::endl;
}

void DSurfTomo:: forward(Eigen::Tensor<float,3> &vs,Eigen::VectorXf &data)
{
    surf.TravelTime(mod,vs,data);
}

void DSurfTomo :: checkerboard()
{
    VectorXf dsyn(surf.num_data);
    forward(vstrue,dsyn);

    // add noise
    for(int i=0;i<surf.num_data;i++){
        //dsyn(i) *= (1.0 + param.noiselevel * gaussian());
        dsyn(i) += param.noiselevel * gaussian();
    }
    surf.obst = dsyn*1.0;
}


void DSurfTomo::inversion(Tensor<float,3> &vsf,VectorXf &dsyn)
{
    int m = surf.num_data;
    int nx = mod.nx, ny = mod.ny;
    int nz = mod.nz;
    int n = (nx-2) * (ny -2) * (nz  - 1);

    // residual vector and velocity variation vector
    VectorXf res(m + n),dv(n);
    res.setZero();
    dv.setZero();

    // compute frechet kernel
    std::string basedir = "kernel";
    VectorXi nonzeros = surf.FrechetKernel(mod,vsf,dsyn,basedir);
    res.segment(0,m) = surf.obst - dsyn;

    // initialize matrix and compute indptr
    int nar = nonzeros.sum();
    csr_matrix<float> smat(m+n,n,nar + n*7);
    smat.indptr[0] = 0;
    for(int i=0;i<m;i++){
        smat.indptr[i+1] = smat.indptr[i] +  nonzeros[i];
    }
    for(int i=m;i<m+n;i++) smat.indptr[i+1] = smat.indptr[m];

    // assembling derivative matrix
    std::cout << "Assembling derivative Matrix ..." << std::endl;
    surf.read_Frechet_Kernel(basedir,smat);
    
    // add regularization terms
    mod.add_regularization(smat,param.smooth);

    // solve equations by lsmr
    std::cout <<"solving linear systems by LSMR ..." << std::endl;
    int itnlim = n * 2; 
    LSMRDict<float> dict(itnlim,10,param.damp,param.smooth,param.nthreads);
    smat.LsmrSolver(res.data(),dv.data(),dict);
    std::cout << "max negative and positive perturbation: " \
                << dv.minCoeff() <<" " << dv.maxCoeff()\
                <<std::endl;
    
    // tackle large variations
    Eigen::TensorMap<Tensor<float,3>>dx(dv.data(),nx-2,ny-2,nz-1);
    float minvel = param.minvel;
    float maxvel = param.maxvel;
    for(int k=0;k<nz-1;k++){
    for(int j=0;j<ny-2;j++){
    for(int i=0;i<nx-2;i++){
        if(dx(i,j,k) > 0.5) dx(i,j,k) = 0.5;
        if(dx(i,j,k) < -0.5) dx(i,j,k) = -0.5;
        float temp = dx(i,j,k) + vsf(i+1,j+1,k);
        if(temp > maxvel) temp = maxvel;
        if(temp < minvel) temp = minvel;
        if(vsf(i+1,j+1,k) !=0.0) vsf(i+1,j+1,k) = temp;
    }}}
}