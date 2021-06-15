#define EIGEN_DONT_PARALLELIZE
#include"SurfaceWave.hpp"
#include"utils.hpp"
#include"surfdisp.hpp"
#include<omp.h>
#include<fstream>
using Eigen::Tensor;
using Eigen::MatrixXd;
using Eigen:: VectorXf;
using Eigen::VectorXd;

// compute period index
int SurfTime :: get_period_index(int period_idx,std::string &wavetype)
{
    int idx;
    if(wavetype == "Rc") idx = period_idx;
    if(wavetype=="Rg") idx = period_idx + kmaxRc;
    if(wavetype=="Lc") idx = period_idx + kmaxRc + kmaxRg;
    if(wavetype == "Lg") idx = period_idx + kmaxRc + kmaxRg + kmaxLc;

    return idx;
}

/**
 * Convert kernels of layer-based model to that of grid-based model
 * @param nz no. of grid points for grid-based model
 * @param rmax no. of layers for layer-based model
 * @param dep,vp,vs,rho grid-based model
 * @param rvp,rvs,rrho rthk layer-based model
 * @param rdcdm parameters sensitivity kernel for layer-based model 
 * @param dcdm parameters sensitivity kernel for grid-based model
 */
void _ParamConvert(float *dep,float *vp,float *vs,float *rho,int nz,
                    int sublayer,float *rthk,float *rvp,float *rvs,
                    float *rrho,int rmax,double *rdcdm,double *dcdm)
{
    // check input parameters
    int nr = nz + (nz-1) * sublayer;
    if(nr !=rmax && vs[0]!= 0.0){
        std::cout << "please check input parameters\n";
        exit(1);
    }

    // convert 
    int n=sublayer + 1;
    if(vs[0]!=0.0)for(int i=0;i<nz;i++){ // no water layers
        double tup = 0.0, tdwn = 0.0;
        int idup = (i-1)*n,iddw = i*n;
        dcdm[i] = 0.0;
        for(int j=0;j<n;j++){
            double tmp = (2.*j+1.0)/(2.*n);
            if(idup +j>=0) tup += rdcdm[idup+j] * tmp;
            if(iddw+j<rmax) tdwn += rdcdm[iddw+j]*(1.0 - tmp);
            if(iddw == rmax-1) tdwn += rdcdm[rmax-1];
        }
        dcdm[i] += tup + tdwn;
    }
    else{ // water layer in the top
        dcdm[0] = rdcdm[0];
        for(int i=1;i<nz;i++){
            dcdm[i] = 0.0;
            double tup = 0.0, tdwn = 0.0;
            int idup = (i-2)*n+1,iddw = (i-1)*n+1;
            for(int j=0;j<n;j++){
                double tmp = (2.*j+1.0)/(2.*n);
                if(idup +j>=1) tup += rdcdm[idup+j] * tmp;
                if(iddw+j<rmax) tdwn += rdcdm[iddw+j]*(1.0 - tmp);
                if(iddw == rmax-1) tdwn += rdcdm[rmax-1];
            }
            dcdm[i] += tup + tdwn;
        }
    }
}

/**
 * Convert kernels of layer-based model to that of grid-based model
 * @param nz no. of grid points for grid-based model
 * @param rmax no. of layers for layer-based model
 * @param dep,vp,vs,rho grid-based model
 * @param rvp,rvs,rrho rthk layer-based model
 * @param rdcdm thk sensitivity kernel for layer-based model 
 * @param dcdm boundary sensitivity kernel for grid-based model
 */
void _DepthConvert(float *dep,float *vp,float *vs,float *rho,int nz,
                    int sublayer,float *rthk,float *rvp,float *rvs,
                    float *rrho,int rmax,double *rdcdh,double *dcdh)
{
    // convert 
    int n=sublayer + 1;
    if(vs[0]!=0.0)for(int i=0;i<nz;i++){ // no water layers
        double tup = 0.0, tdwn = 0.0;
        int idup = (i-1)*n,iddw = i*n;
        dcdh[i] = 0.0;
        for(int j=0;j<n;j++){
            double tmp = 1.0 / n;
            if(idup +j>=0) tup += rdcdh[idup+j] * tmp;
            if(iddw+j<rmax) tdwn += -rdcdh[iddw+j]*tmp;
        }
        dcdh[i] += tup + tdwn;
    }
    else{ // water layer in the top
        dcdh[0] = -rdcdh[0];
        for(int i=1;i<nz;i++){
            double tup = 0.0, tdwn = 0.0;
            int idup = (i-2)*n+1,iddw = (i-1)*n+1;
            dcdh[i] = 0.0;
            for(int j=0;j<n;j++){
                double tmp = 1.0 / n;
                if(idup +j>=1) tup += rdcdh[idup+j] * tmp;
                if(iddw+j<rmax) tdwn += -rdcdh[iddw+j]*tmp;
                if(idup + j ==0) tup += rdcdh[0];
            }
            dcdh[i] += tup + tdwn;
        }
    }
}

/**
 * refine grid based model to layerd based model, water layers considered
 * 
 * \param dep,vp,vs,rho grid-based model
 * \param nz no. of points in grid-based model
 * \param sublayer insert sublayer points to construct layered model
 * \param rdep, rvp, rvs, rrho, rthk layered velocity model 
 * 
 * \return rmax no. of layers of refined model
 */
int grid2LayerModel(float *dep,float *vp,float *vs,float *rho,int nz,
                    int sublayers,float *rthk,float *rvp,float *rvs,
                    float *rrho)
{
    float thk;
    int n = sublayers + 1;
    int rmax,istart;
    if(vs[0] == 0.0){ // tackle water layer
        rmax = nz-1 + (nz-2) * sublayers + 1;
        istart = 1;
        rthk[0] = dep[1] - dep[0];
        rrho[0] = rho[0];
        rvs[0] = 0.0;
        rvp[0] = vp[0];
    }
    else{
        rmax = nz + (nz-1) * sublayers;
        istart = 0;
    }

    int k = istart;
    for(int i=istart;i<nz-1;i++){
        thk = (dep[i+1] - dep[i]) / n;
        // parameters in each layer are midpoint values in grid-based model
        for(int j=0;j<n;j++){ 
            rthk[k] = thk;
            rvp[k] = vp[i] + (2*j+1)*(vp[i+1]-vp[i])/(2.0*n);
            rvs[k] = vs[i] + (2*j+1)*(vs[i+1]-vs[i])/(2.0*n);
            rrho[k] = rho[i] + (2*j+1)*(rho[i+1]-rho[i])/(2.0*n);
            k++ ;
        }
    }

    // half space
    k = rmax-1;
    rthk[k] = 0.0;
    rvp[k] = vp[nz-1];
    rvs[k] = vs[nz-1];
    rrho[k] = rho[nz-1];

    return rmax;
}

/**
 * internal function to compute disperion for one cell
 * Parameters:
 * ====================================================
 *  vs,vp,rho,dep : parameters for one cell
 *  nz            : no. of points in vertical direction
 *  sublayer      : insert sublayer layers in adjacent points
 *  tRc...        : period vector for each wavetype
 *  vmap :dispersion for one cell, size(kmaxRc+kmaxRg+kmaxLc+kmaxLg)
*/
void _disp1D(float *vs,float *vp,float *rho,float *dep,int nz,
            int sublayer,VectorXd &tRc,VectorXd &tRg,
            VectorXd &tLc,VectorXd &tLg,double *vmap)
{
    // convert grid-based model to layer-based model
    int rmax = nz + (nz-1) * sublayer;
    float *rvp = new float [rmax],*rvs = new float [rmax];
    float *rthk = new float [rmax],*rrho = new float [rmax];
    rmax = grid2LayerModel(dep,vp,vs,rho,nz,sublayer,rthk,rvp,rvs,rrho);

    // prepare parameters
    int kmaxRc = tRc.size(),kmaxRg = tRg.size();
    int kmaxLc= tLc.size(), kmaxLg = tLg.size();
    bool sphere = true, keep_flat = false;
    int mode = 0;

    // compute dispersion
    if(kmaxRc > 0){
        surfdisp(rthk,rvp,rvs,rrho,rmax,tRc.data(),vmap,
                kmaxRc,"Rc",mode,sphere,keep_flat);
    }
    if(kmaxRg > 0){
        _RayleighGroup(rthk,rvp,rvs,rrho,rmax,tRg.data(),
                        vmap + kmaxRc,kmaxRg,mode,sphere);
    }
    if(kmaxLc > 0){
        surfdisp(rthk,rvp,rvs,rrho,rmax,tLc.data(),vmap + kmaxRc+kmaxRg,
                 kmaxLc,"Lc",mode,sphere,keep_flat);
    }
    if(kmaxLg > 0){
        _LoveGroup(rthk,rvs,rrho,rmax,tLg.data(),vmap+kmaxRc+kmaxRg+kmaxLc,
                    kmaxLg,mode,sphere);
    }

    // free space
    delete[] rrho;
    delete[] rvp;
    delete[] rvs;
    delete[] rthk;
}

/**
 * internal function to compute dispersion and sensitivity kernels
 * Parameters:
 * ====================================================
 *  vs,vp,rho,dep : parameters for one cell
 *  nz            : no. of points in vertical direction
 *  sublayer      : insert sublayer layers in adjacent points
 *  tRc...        : period vector for each wavetype
 *  vmap          : dispersion for one cell, size(kmaxRc+kmaxRg+kmaxLc+kmaxLg)
 *  kvs/vp/rho    : sensitivity kernel for vs/vp/rho, shape(nz,nt),column major
*/
void _Kernel1D(float *vs,float *vp,float *rho,float *dep,int nz,
            int sublayer,VectorXd &tRc,VectorXd &tRg,VectorXd &tLc,
            VectorXd &tLg,double *vmap,double *kvs,double *kvp,
            double *krho)
{
    // convert grid-based model to layer-based model
    int rmax = nz + (nz-1) * sublayer;
    float *rvp = new float [rmax],*rvs = new float [rmax];
    float *rthk = new float [rmax],*rrho = new float [rmax];
    rmax = grid2LayerModel(dep,vp,vs,rho,nz,sublayer,rthk,rvp,rvs,rrho);

    // allocate space for kernel
    double *ur = new double [rmax], *uz = new double[rmax]; 
    double *tr = new double [rmax], *tz = new double[rmax];
    double *dudar = new double[rmax], *dudbr = new double[rmax];
    double *dudrr = new double [rmax],*dudhr = new double [rmax];
    double *dcdar = new double[rmax], *dcdbr = new double[rmax];
    double *dcdrr = new double [rmax],*dcdhr = new double [rmax];

    // prepare parameters
    int kmaxRc = tRc.size(),kmaxRg = tRg.size();
    int kmaxLc= tLc.size(), kmaxLg = tLg.size();
    bool sphere = true;
    int mode = 0, iflsph = 1;

    // compute dispersion and sensitivity kernel
    if(kmaxRc > 0){ // rayleigh group
        surfdisp(rthk,rvp,rvs,rrho,rmax,tRc.data(),vmap,
                kmaxRc,"Rc",mode,sphere);
        for(int i=0;i<kmaxRc;i++){
            double cg;
            sregn96_(rthk,rvp,rvs,rrho,rmax,tRc.data()+i,vmap+i,&cg,
                    ur,uz,tr,tz,dcdar,dcdbr,dcdhr,dcdrr,iflsph);
            _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dcdar,kvp+i*nz);
            _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dcdbr,kvs+i*nz);
            _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dcdar,krho+i*nz);          
        }   
    }
    if(kmaxRg > 0){
        double cp[kmaxRg],t1[kmaxRg],t2[kmaxRg],c1[kmaxRg],c2[kmaxRg];
        double dt = 0.01;
        for(int i=0;i<kmaxRg;i++){
            t1[i] = tRg(i) * (1.0 + 0.5 * dt);
            t2[i] = tRg(i) * (1.0 - 0.5 * dt);
        }
        surfdisp(rthk,rvp,rvs,rrho,rmax,tRg.data(),cp,kmaxRg,"Rc",mode,sphere);
        surfdisp(rthk,rvp,rvs,rrho,rmax,t1,c1,kmaxRg,"Rc",mode,sphere);
        surfdisp(rthk,rvp,rvs,rrho,rmax,t2,c2,kmaxRg,"Rc",mode,sphere);
        for(int i=0;i<kmaxRg;i++){
            int k = i + kmaxRc;
            sregnpu_(rthk,rvp,rvs,rrho,rmax,tRg.data()+i,cp+i,vmap+k,
                     ur,uz,tr,tz,t1+i,c1+i,t2+i,c2+i,dcdar,
                     dcdbr,dcdhr,dcdrr,dudar,dudbr,dudhr,dudrr,iflsph);
            _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dudar,kvp+k*nz);
            _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dudbr,kvs+k*nz);
            _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dudrr,krho+k*nz);          
        } 
    }
    if(kmaxLc > 0){
        surfdisp(rthk,rvp,rvs,rrho,rmax,tLc.data(),vmap+kmaxRc+kmaxRg,
                kmaxLc,"Lc",mode,sphere);
        for(int i=0;i<kmaxLc;i++){
            double cg;
            int k = i + kmaxRc + kmaxRg;
            slegn96_(rthk,rvs,rrho,rmax,tLc.data(),vmap+k,
                    &cg,ur,tr,dcdbr,dcdhr,dcdrr,iflsph);
            _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dcdbr,kvs+k*nz);
            _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dcdrr,krho+k*nz);          
        }  
    }
    if(kmaxLg > 0){
        double cp[kmaxLg],t1[kmaxLg],t2[kmaxLg],c1[kmaxLg],c2[kmaxLg];
        double dt = 0.01;
        for(int i=0;i<kmaxLg;i++){
            t1[i] = tLg(i) * (1.0 + 0.5 * dt);
            t2[i] = tLg(i) * (1.0 - 0.5 * dt);
        }
        surfdisp(rthk,rvp,rvs,rrho,rmax,tLg.data(),cp,
                kmaxLg,"Lc",mode,sphere);
        surfdisp(rthk,rvp,rvs,rrho,rmax,t1,c1,
                kmaxLg,"Lc",mode,sphere);
        surfdisp(rthk,rvp,rvs,rrho,rmax,t2,c2,
                kmaxLg,"Lc",mode,sphere);
        for(int i=0;i<kmaxLg;i++){
            int k = i + kmaxRc + kmaxRg + kmaxLc;
            slegnpu_(rthk,rvs,rrho,rmax,tLg.data()+i,cp+i,vmap+k,
                    ur,tr,t1+i,c1+i,t2+i,c2+i,dcdbr,dcdhr,dcdrr,
                    dudbr,dudhr,dudrr,iflsph);
            _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dudbr,kvs+k*nz);
            _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dudrr,krho+k*nz);          
        } 
    }

    // free space
    delete[] rrho;delete[] rvp;delete[] rvs;delete[] rthk;
    delete[] ur; delete [] uz; delete [] tr; delete[] tz;
    delete [] dudrr; delete[] dudar; delete[] dudbr; delete[] dudhr;
    delete [] dcdrr; delete[] dcdar; delete[] dcdbr; delete[] dcdhr;
}

/**
 * Compute dispersion map
 * Parameters:
 * ===========================================
 *      mod   : background model
 *      vs    : vs velocity (nx,ny,nz)
 *      vmap  : dispersion map(nx*ny,nt), 
 *              fundamental mode Rc,Rg,Lc,Lg, and first mode ...
*/
void SurfTime ::DipersionMap(MOD3d &mod,Eigen::Tensor<float,3> &vs,Eigen::MatrixXd &vmap)
{
    // get vs-dimension
    int nx = mod.nx,ny = mod.ny,nz = mod.nz;
    int nt = vmap.cols();

    // copy vs to vs0
    MatrixXd vmap0(nt,nx*ny);
    Tensor<float,3> vs0(nz,nx,ny);
    for(int i=0;i<ny;i++){
    for(int j=0;j<nx;j++){
    for(int k=0;k<nz;k++){
        vs0(k,j,i) = vs(j,i,k);
    }}}

    // compute dispersion map
    int nthreads = this->num_threads;
    omp_set_num_threads(nthreads);
    #pragma omp parallel for shared(vs0,vmap0,mod)
    for(int n=0;n<nx*ny;n++){
        int i = n % nx;
        int j = n/nx;

        // convert vs to vp, rho by empirical relations
        float v[nz],rho[nz],vp[nz];
        for(int k=0;k<nz;k++){
            v[k] = vs0(k,i,j);
            mod.empirical_relation(v[k],vp[k],rho[k]);
        }

        // compute frechet kernel
        _disp1D(v,vp,rho,mod.dep.data(),nz,sublayer,
                tRc,tRg,tLc,tLg,&vmap0(0,n));
    }

    vmap = vmap0.transpose();
}

/**
 * Compute 1-D surface wave frechet kernel
 * Parameters:
 *  mod     : 3-D initial model
 *  vs      : 3-D shear wave, shape(nx,ny,nz)
 *  vmap    : dispersion map, shape(nx*ny,nt)
 *  kernel  : frechet kernel, shape(nx*ny,nz,nt) 
 */
void SurfTime:: SurfWaveKernel(MOD3d &mod,Tensor<float,3> &vs,
                                MatrixXd &vmap,Tensor<double,3> &kernel)
{
    // get vs-dimension
    int nx = mod.nx,ny = mod.ny,nz = mod.nz;
    int nt = vmap.cols();

    // copy vs to vs0 
    Tensor<float,3> vs0(nz,nx,ny);
    for(int i=0;i<ny;i++){
    for(int j=0;j<nx;j++){
    for(int k=0;k<nz;k++){
        vs0(k,j,i) = vs(j,i,k);
    }}}

    // compute frechet kernel and dispersion map
    MatrixXd vmap0(nt,nx*ny);
    Tensor<double,3> kernel0(nz,nt,nx*ny);
    int nthreads = this->num_threads;
    omp_set_num_threads(nthreads);
    #pragma omp parallel for shared(vs0,vmap0,mod,kernel0)
    for(int n=0;n<nx*ny;n++){
        int i = n % nx;
        int j = n/nx;

        // convert vs to vp, rho by empirical relations
        float v[nz],rho[nz],vp[nz];
        for(int k=0;k<nz;k++){
            v[k] = vs0(k,i,j);
            mod.empirical_relation(v[k],vp[k],rho[k]);
        }

        // compute dispersion map for vs,vp and rho
        Eigen::MatrixXd kvs(nz,nt),kvp(nz,nt),krho(nz,nt);
        kvs.setZero(); kvp.setZero();krho.setZero();
        _Kernel1D(v,vp,rho,mod.dep.data(),nz,sublayer,
                tRc,tRg,tLc,tLg,&vmap0(0,n),kvs.data(),
                kvp.data(),krho.data());

        // convert to vs by empirical relations
        for(int ii=0;ii<nt;ii++){
            float drda,dadb;
        for(int jj=0;jj<nz;jj++){
            if(v[jj] == 0.0 && jj == 0) continue; // water layer
            mod.empirical_deriv(vp[jj],v[jj],drda,dadb);
            kernel0(jj,ii,n) = kvs(jj,ii) + kvp(jj,ii) * dadb 
                              + krho(jj,ii) * drda * dadb;
        }}
    }

    // copy pv0 and kernel0 to pv and kernel
    vmap = vmap0.transpose();
    for(int i=0;i<nt;i++){
    for(int j=0;j<nz;j++){
    for(int k=0;k<nx*ny;k++){
        kernel(k,j,i) = kernel0(j,i,k);
    }}}
    
}

void SurfTime :: TravelTime(MOD3d &mod,Eigen::Tensor<float,3> &vs,Eigen::VectorXf &data)
{
    // get vs-dimension
    int nx = mod.nx,ny = mod.ny;
    MatrixXd pv(nx*ny,kmax);

    // compute dispersion map
    DipersionMap(mod,vs,pv);

    // get traveltime for every source-receiver pair
    int np = Pairs.size();
    int nthreads = this->num_threads;
    omp_set_num_threads(nthreads);
    #pragma omp parallel for shared(vs,pv,mod,data)
    for(int n=0;n<np;n++){
        StationPair &p = Pairs[n];
        int nr = p.nr;
        int idx = get_period_index(p.period_idx,p.wavetype);
        float ttime[nr];

        synthetic(nx,ny,&pv(0,idx),mod.goxd,mod.gozd,mod.dvxd,mod.dvzd,p.srcx,p.srcz,
                &p.rcx[0],&p.rcz[0],nr,ttime);

        for(int i=0;i<nr;i++){
            data(p.counter + i) = ttime[i];
        }
    }
}

Eigen::VectorXi 
SurfTime :: FrechetKernel(MOD3d &mod,Tensor<float,3> &vs,VectorXf &data,
                            std::string &save_dir)
{
  // get vs-dimension
    int nx = mod.nx,ny = mod.ny,nz= mod.nz;
    MatrixXd pv(nx*ny,kmax);
    Tensor<double,3> kernel(nx*ny,nz,kmax);

    // compute dispersion map and frechet kernel
    std::cout << "computing Surface Wave 1-D Frechet Kernel ..." << std::endl;
    SurfWaveKernel(mod,vs,pv,kernel);

    // --------------------------------------------------
    // get traveltime for every source-receiver pair
    // ---------------------------------------------------
    std::cout << "computing traveltimes and 2-D Frechet Kernel ..." <<std::endl;
    int np = Pairs.size(); // no. of pairs
    
    // compute no. of non-zero elements in frechet matrix for each data
    Eigen::VectorXi nonzeros(num_data);

    //std::vector<int> nonzeros(nthreads);
    std::string cmd = "mkdir -p " + save_dir;
    int flag= system(cmd.c_str()); // mkdir
    
    // parallelization
    int nthreads = this->num_threads;
    omp_set_num_threads(nthreads);
    #pragma omp parallel for shared(vs,pv,mod,data,kernel)
    for(int nn=0;nn<nthreads;nn++){
        // open a file to save frechet kernel
        std::string filename = save_dir + "/" + std::to_string(nn) +".txt";
        FILE *fp;
        if((fp=fopen(filename.c_str(),"w"))==NULL){
            std::cout << "cannot open file "<< filename << std::endl;
            exit(0);
        }

        // loop around all station pairs
        for(int pp=nn;pp<np;pp+=nthreads){
            StationPair &p = Pairs[pp];
            int nr = p.nr;
            int idx = get_period_index(p.period_idx,p.wavetype);
            float ttime[nr];

            // compute 3-D pseudo sensitivity kernel
            int n = (nx-2)*(ny-2)*(nz-1);
            Eigen::MatrixXf frechet(n,nr);
            CalFrechet(nx,ny,nz,&pv(0,idx),mod.goxd,mod.gozd,
                        mod.dvxd,mod.dvzd,p.srcx,p.srcz,&p.rcx[0],
                        &p.rcz[0],nr,ttime,&kernel(0,0,idx),frechet.data());
            
            // loop around all receivers
            int counter = p.counter;
            for(int i=0;i<nr;i++){
                int nar = (frechet.col(i).array().abs()>0.0).cast<int>().sum();
                nonzeros(counter) = nar;
                fprintf(fp,"# %d %d\n",counter,nar);
                for(int j=0;j<n;j++){
                    float tmp = frechet(j,i);
                    if(std::abs(tmp)> 0.0) {
                        fprintf(fp,"%d %g\n",j,tmp);
                    }
                }
                data(p.counter + i) = ttime[i];
                counter += 1;
            }
        }
        fclose(fp);
    }
    return nonzeros;
}

void SurfTime::
read_Frechet_Kernel(std::string &basedir,csr_matrix<float> &smat)
{
    int nthreads = this->num_threads;
    const int MAX_LEN = 100;
    omp_set_num_threads(nthreads);
    #pragma omp parallel for shared(smat)
    for(int p=0;p<nthreads;p++){
        std::string filename = basedir + "/" +  std::to_string(p) + ".txt";

        // open file
        char line[MAX_LEN];
        FILE *fp;
        if((fp=fopen(filename.c_str(),"r"))==NULL){
            std::cout << "cannot open file "<< filename << std::endl;
            exit(0);
        }

        int rw_idx,nar;
        char dummy;
        while(fgets(line,sizeof(line),fp)!=NULL){
            if(line[0] != '#') break;
            sscanf(line,"%c%d%d",&dummy,&rw_idx,&nar);
            int start = smat.indptr[rw_idx];
            int end = smat.indptr[rw_idx + 1];

            if(end - start != nar){
                std::cout << "format error!\n";
                exit(1);
            }

            for(int i=start;i<end;i++){
                assert(fscanf(fp,"%d%f\n",smat.indices + i,smat.data + i) == 2);
            }
        }

        //close and delete file 
        fclose(fp);
        int ierr = system(("rm -r "+ filename).c_str() );
    } 
}

void SurfTime:: 
write_disper(Eigen::VectorXf &dsyn,std::string &filename){
    std::ofstream outfile;
    outfile.open(filename);

    // output
    int c = 0; // counter
    const double rad2deg = 180. / M_PI;
    for(unsigned int n=0;n<Pairs.size();n++){
        StationPair &p = Pairs[n];
        float elat = (M_PI * 0.5 - p.srcx) * rad2deg, elon = p.srcz * rad2deg;
        float T;
        if(p.wavetype == "Rc")
            T = tRc(p.period_idx);
        else if(p.wavetype == "Rg")
            T = tRg(p.period_idx);
        else if(p.wavetype == "Lc")
            T = tLc(p.period_idx);
        else 
            T = tLg(p.period_idx);
        for(int i=0;i<p.nr;i++){
            outfile << sta_dist[c] << " " << obst(c) << " " << dsyn(c) << " ";
            outfile << elat << " " << elon << " ";
            float  slat = (M_PI*0.5 - p.rcx[i]) * rad2deg, slon = p.rcz[i] * rad2deg;
            outfile << slat << " " << slon << " " << T << " " << p.wavetype << "\n";
            c += 1;
        }
    }
    outfile.close();
}