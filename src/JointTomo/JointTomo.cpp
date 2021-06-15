#define EIGEN_DONT_PARALLELIZE
#include"IOFunction.hpp"
#include"tomography.hpp"
#include"utils.hpp"
#include<omp.h>
using Eigen::Tensor;
using Eigen::VectorXf;
using Eigen::VectorXi;
const double DEG2RAD = M_PI/180.;

/**
 * read all input parameters, observations and initial model
 * @param paramfile file contains parameter
 * @param surfdata file contains observed dispersion data
 * @param gravdata file contains gravity data
 * @param modfile initial model
 * @param modtrue True model, for synthetic test 
 */
void JointTomo::
readdata(std::string &paramfile,std::string &modfile,std::string &surfdata,
        std::string &gravdata,std::string &gravmat,std::string &refmod,
        std::string &modtrue)
{
   // step1 : read paramfile
    std::ifstream infile;
    std::string line;
    std ::istringstream info;    
    infile.open(paramfile);

    // read parfile
    param.nthreads = __read_parfile(infile,paramfile,mod,surf,param);

    // read synthetic test params
    int remove_average;
    skipread(infile,line,"%f%f",&param.noiselevel,&param.noiselevel1);
    skipread(infile,line,"%f",&param.p);
    skipread(infile,line,"%f%f",&param.weight1,&param.weight2);
    skipread(infile,line,"%d",&remove_average);
    param.remove_average = remove_average == 1;
    infile.close();

    // allocate data parameters
    int nx = mod.nx,ny = mod.ny, nz= mod.nz;
    this->unknowns = (nx-2) * (ny -2) * (nz -1 );
    if(param.ifsyn) vstrue.resize(nx,ny,nz);

    // step2 :read surface wave traveltime data
    __read_traveltime(surfdata,surf);

    // step3 : read gravity data and remove average value
    printf("\nGravity Data:\n");
    printf("===================================\n");
    obsg.read_obs_data(gravdata);
    std::cout <<"The number of gravity measurements = " <<obsg.np << std::endl;
    std::cout << "Remove average value of synthetic data = ";
    if(param.remove_average){
        std::cout << "True \n";
    }
    else{
        std::cout << "False \n";
    }
    this->num_data = surf.num_data + obsg.np;

    // step4: read initial model
    __read_InputModel(modfile,mod.vs,mod.dep.data(),true);

    // step5 : read ture model if required
    if(param.ifsyn == 1){
        float z[nz];
        __read_InputModel(modtrue,vstrue,z);
    }

    // step6: read refmodel
    modref = mod;
   if(refmod!="None"){
        std::cout<<"\nThe gravity reference model is " + refmod <<std::endl;
        float z[nz];
        __read_InputModel(refmod,modref.vs,z);
    }
    else{
        std::cout<<"\nNo reference model is given, so the average \
                    of initial model is used" << std::endl;
        for(int k=0;k<nz;k++){
            float mean = 0.0;
            for(int j=0;j<ny-2;j++){
            for(int i=0;i<nx-2;i++){
                mean += mod.vs(i+1,j+1,k);
            }}
            mean /= (nx-2) * (ny - 2);
            for(int j=0;j<ny;j++){
            for(int i=0;i<nx;i++){
                modref.vs(i,j,k) = mean;
            }}
        }
    }
    
    //step7 : read gravity matrix
    std::cout << "Reading gravity matrix ...\n" << std::endl;
    int nar = 0;
    FILE *pin;
    if((pin=popen(("grep -v '#' " + gravmat + "| wc -l").c_str(), "r"))==NULL){
        std::cout << "cannot open file "<< gravmat << std::endl;
        exit(0);
    }
    assert(fscanf(pin,"%d",&nar)==1);
    pclose(pin);
    gmat.initialize(obsg.np,unknowns,nar);
    gmat.read(gravmat);  
} 

void JointTomo ::forward(Tensor<float,3> &vs,VectorXf &dsyn,VectorXf &dg)
{      
    surf.TravelTime(mod,vs,dsyn);
    modref.gravity(gmat,vs,dg);
} 
  
void JointTomo:: checkerboard()
{
    VectorXf dsyn(surf.num_data),dg(obsg.np);
    forward(vstrue,dsyn,dg);

    // add noise
    for(int i=0;i<surf.num_data;i++){
        //surf.obst(i) = dsyn(i) * (1.0 +  param.noiselevel * gaussian());
        surf.obst(i) = dsyn(i) + param.noiselevel * gaussian();
    }

    for(int i=0;i<obsg.np;i++){
        //obsg.Gr[i] = dg(i) * (1.0 + param.noiselevel1 * gaussian());
        obsg.Gr[i] = dg(i) + param.noiselevel1 * gaussian();
    }
}

/**
 * Assemble gravity and SWD matrix to a global one, and add weights
 * -------------------------------------------------------------------
 * empirical relation is used to convert dg/drho to dg/dvs
 * convert change of density to change of S wave velocity
 * some other empirical relation maybe better
 * here we use:
 * rho = 1.6612*vp-0.4721*vp**2 + 
 *          0.0671*vp**3 - 0.0043*vp**4 + 0.000106*vp**5
 * vp = 0.9409 + 2.0947*vs - 0.8206*vsz**2+ 0.2683*vs**3 - 0.0251*vs**4
 * then G*drho = g0 could be changed to G'*dvs = g0 
 */ 
void JointTomo:: 
assemble(Tensor<float,3> &vsf,csr_matrix<float> &smat,
         float weight1,float weight2)
{
    // get model and matrix dimensions
    int nx = mod.nx, ny = mod.ny;
    int ngrav = gmat.rows();
    int nsurf = surf.num_data;

    // assemble
    for(int r = 0;r<ngrav;r++){ // loop around all rows
        int start = gmat.indptr[r];
        int end = gmat.indptr[r + 1];
        int rwc = nsurf + r; 
        smat.indptr[rwc + 1] = smat.indptr[rwc] + end - start;
    }
    int nthreads = surf.num_threads;  
    omp_set_num_threads(nthreads);
    #pragma omp parallel for shared(smat,vsf)
    for(int r=0;r<ngrav;r++){
        int start = gmat.indptr[r];
        int end = gmat.indptr[r + 1];
        int rwc = nsurf + r;
        for(int c=start;c<end;c++){ // loop around nonzero columns
            int col = gmat.indices[c];
            int k = col /((nx-2) * (ny-2));
            int ridx = col %((nx-2) * (ny-2));
            int i = ridx%(nx-2),j = ridx/(nx-2);

            // compute drho/dvs, store it in tmp1 * tmp2
            float vp,vs,rho;
            vs = vsf(i,j,k);
            mod.empirical_relation(vs,vp,rho);
            float tmp1,tmp2;
            mod.empirical_deriv(vp,vs,tmp1,tmp2);

            // append gmat to smat
            int count = smat.indptr[rwc] + c - start;
            smat.data[count] = gmat.data[c] * tmp1 * tmp2;
            smat.indices[count] = gmat.indices[c];
        }
    }

    // add weights
    for(int i=0;i<nsurf+ngrav;i++){
    for(int j=smat.indptr[i];j<smat.indptr[i+1];j++){
        if(i < nsurf) 
            smat.data[j] *= weight1;
        else 
            smat.data[j] *= weight2;
    }}
}

void JointTomo:: inversion(Tensor<float,3> &vsf,VectorXf &dsyn,VectorXf &dg)
{
    int nx = mod.nx, ny = mod.ny;
    int nz = mod.nz;
    int m = num_data,n = unknowns;
    VectorXf res(m + n);

    // compute frechet kernel, and gravity anomaly
    std::string basedir = "kernelJ";
    VectorXi nonzeros = surf.FrechetKernel(mod,vsf,dsyn,basedir);
    modref.gravity(gmat,vsf,dg);
    if(param.remove_average){ // remove average value if needed
        float mean = dg.sum() / dg.size();
        dg.array() -= mean;
    }

    // renew residuals vector
    res.setZero();
    Eigen::Map<VectorXf> Gr(obsg.Gr,obsg.np);
    res.segment(0,surf.num_data) = surf.obst - dsyn;
    res.segment(surf.num_data,obsg.np) = Gr - dg;

    // compute weights according to Julia(2000)
    float sigma1,sigma2;
    if(param.ifsyn == 1){
        //sigma1 = mean * (1. + param.noiselevel) / sqrt(surf.num_data);
        //sigma2= meang * (1. + param.noiselevel1) / sqrt(surf.num_data);
        sigma1 = param.noiselevel / sqrt(surf.num_data);
        sigma2 = param.noiselevel1 / sqrt(surf.num_data);
    }
    else {
        sigma1 = param.weight1 / sqrt(surf.num_data);
        sigma2 = param.weight2 / sqrt(surf.num_data);
    }
    
    sigma1 = sqrt(param.p /surf.num_data) / sigma1;
    sigma2 = sqrt((1-param.p)/obsg.np) / sigma2;

    // initialize global sparse matrix  
    int nar = nonzeros.sum(); 
    csr_matrix<float> smat(m+n,n,nar + gmat.nonzeros + n * 7);
    smat.indptr[0] = 0;
    for(int i=0;i<surf.num_data;i++){
        smat.indptr[i+1] = smat.indptr[i] +  nonzeros[i];
    }
    for(int i=surf.num_data;i<m+n;i++) smat.indptr[i+1] = smat.indptr[i];

    // assembling derivative matrix
    std::cout << "Assembling derivative Matrix ..." << std::endl;
    surf.read_Frechet_Kernel(basedir,smat);
    assemble(vsf,smat,sigma1,sigma2);

    // add weights
    for(int i=0;i<surf.num_data;i++){
        res(i) *= sigma1;
    }
    for(int i=0;i<obsg.np;i++){
        res(i+surf.num_data) *= sigma2;
    }

    // add regularization terms
    mod.add_regularization(smat,param.smooth);

    // solve equations by lsmr
    std::cout <<"solving linear systems by LSMR ..." << std::endl;
    int itnlim = n * 2;
    VectorXf dv(n);
    LSMRDict<float> dict(itnlim,10,param.damp,param.smooth,param.nthreads);
    smat.LsmrSolver(res.data(),dv.data(),dict);
    std::cout << "max negative and positive perturbation: " \
                << dv.minCoeff() <<" " << dv.maxCoeff()\
                <<std::endl;

    // renew vsf and tackle large variations
    Eigen::TensorMap<Tensor<float,3>>dx(dv.data(),nx-2,ny-2,nz-1);
    float minvel = param.minvel;
    float maxvel = param.maxvel;
    for(int k=0;k<nz-1;k++){
    for(int j=0;j<ny-2;j++){
    for(int i=0;i<nx-2;i++){
    	float temp = dx(i,j,k);
        if(temp > 0.5) temp = 0.5;
        if(temp < -0.5) temp = -0.5;
        dx(i,j,k) = temp;

        temp = temp + vsf(i+1,j+1,k);
        if(temp > maxvel) temp = maxvel;
        if(temp < minvel) temp = minvel;
        if(vsf(i+1,j+1,k) !=0.0) vsf(i+1,j+1,k) = temp;
    }}}
}
