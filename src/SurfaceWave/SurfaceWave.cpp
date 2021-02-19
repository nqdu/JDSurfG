#define EIGEN_DONT_PARALLELIZE
#include"SurfaceWave.hpp"
#include"utils.hpp"
#include"surfdisp.hpp"
#include"openmp.hpp"
#include<fstream>
using Eigen::Tensor;
using Eigen::MatrixXd;
using Eigen:: VectorXf;
using Eigen::VectorXd;

// compute period index
int SurfTime :: get_period_index(int period_idx,std::string wavetype)
{
    int idx;
    if(wavetype == "Rc") idx = period_idx;
    if(wavetype=="Rg") idx = period_idx + kmaxRc;
    if(wavetype=="Lc") idx = period_idx + kmaxRc + kmaxRg;
    if(wavetype == "Lg") idx = period_idx + kmaxRc + kmaxRg + kmaxLc;

    return idx;
}

/**
 * refine grid based model to layerd based model
 * 
 * Parameters:
 * -------------------------------------------------
 *  sublayer: no. of grids inserted in each raw layer
 *  mmax: number of depth grid points in the model
 *  dep, vp, vs, rho: the depth-grid model parameters
 *  rmax: number of layers in the fined layered model
 * 
 * Returns:
 * --------------------------------------------------
 *  rdep, rvp, rvs, rrho, rthk: the refined layered velocity model 
 */
void grid2LayerModel(float *dep,float *vp,float *vs,float *rho,int nz,
                    int sublayers,float *rthk,float *rvp,float *rvs,float *rrho)
{
    float thk;
    int n = sublayers + 1;
    int rmax = nz + (nz-1) * sublayers;

    for(int i=0;i<nz-1;i++){
        thk = (dep[i+1] - dep[i]) / n;
        for(int j=0;j<n;j++){
            int k = i*n + j;
            rthk[k] = thk;
            rvp[k] = vp[i] + (2*j+1)*(vp[i+1]-vp[i])/(2.0*n);
            rvs[k] = vs[i] + (2*j+1)*(vs[i+1]-vs[i])/(2.0*n);
            rrho[k] = rho[i] + (2*j+1)*(rho[i+1]-rho[i])/(2.0*n);
        }
    }

    // half space
    int k = rmax-1;
    rthk[k] = 0.0;
    rvp[k] = vp[nz-1];
    rvs[k] = vs[nz-1];
    rrho[k] = rho[nz-1];
}

/**
 * internal function to compute disperion for one cell
 * Parameters:
 *  vmap    : dispersion for one cell, size(kmaxRc+kmaxRg+kmaxLc+kmaxLg)
 *  vs,vp,rho,dep : parameters for one cell
 *  nz            : no. of points in vertical direction
 *  sublayer      : insert sublayer layers in adjacent points
 *  tRc...        : period vector for each wavetype
*/
void __dispMap__(double *vmap,float *vs,float *vp,float *rho,float *dep,int nz,
            int sublayer,VectorXd &tRc,VectorXd &tRg,VectorXd &tLc,VectorXd &tLg)
{
    // convert grid-based model to layer-based model
    int rmax = nz + (nz-1) * sublayer;
    float *rvp = new float [rmax],*rvs = new float [rmax];
    float *rthk = new float [rmax],*rrho = new float [rmax];
    grid2LayerModel(dep,vp,vs,rho,nz,sublayer,rthk,rvp,rvs,rrho);

    // map vmap to eigen's vector
    int kmaxRc = tRc.size(),kmaxRg = tRg.size();
    int kmaxLc= tLc.size(), kmaxLg = tLg.size();
    Eigen::Map<VectorXd> pv(vmap,kmaxLc+kmaxLg+kmaxRc+kmaxRg);

    // compute dispersion
    if(kmaxRc > 0){
        Eigen::VectorXd pvRc(kmaxRc);
        pvRc.setZero();
        surfdisp(rthk,rvp,rvs,rrho,rmax,tRc.data(),pvRc.data(),kmaxRc,"Rc");
        pv.segment(0,kmaxRc) = pvRc;
    }
    if(kmaxRg > 0){
        Eigen::VectorXd pvRc(kmaxRg);
        surfdisp(rthk,rvp,rvs,rrho,rmax,tRg.data(),pvRc.data(),kmaxRg,"Rg");
        pv.segment(kmaxRc,kmaxRg) = pvRc;
    }
    if(kmaxLc > 0){
        Eigen::VectorXd pvRc(kmaxLc);
        surfdisp(rthk,rvp,rvs,rrho,rmax,tLc.data(),pvRc.data(),kmaxLc,"Lc");
        pv.segment(kmaxRc+kmaxRg,kmaxLc) = pvRc;
    }
    if(kmaxLg > 0){
        Eigen::VectorXd pvRc(kmaxLg);
        surfdisp(rthk,rvp,rvs,rrho,rmax,tLg.data(),pvRc.data(),kmaxLg,"Lg");
        pv.segment(kmaxRc+kmaxRg+kmaxLc,kmaxLg) = pvRc;
    }

    // free space
    delete[] rrho;
    delete[] rvp;
    delete[] rvs;
    delete[] rthk;
}

/**
 * Compute dispersion map
 * Parameters:
 *      mod         : background model
 *      vs          : vs velocity (nx,ny,nz)
 *      pv          : dispersion map(nx*ny,nt), 
 *                    fundamental mode Rc,Rg,Lc,Lg, and first mode ...
*/
void SurfTime ::DipersionMap(MOD3d &mod,Eigen::Tensor<float,3> &vs,Eigen::MatrixXd &pv)
{
    // get vs-dimension
    int nx = mod.nx,ny = mod.ny,nz = mod.nz;
    int nt = pv.cols();

    // copy vs to vs0
    MatrixXd pv0(nt,nx*ny);
    Tensor<float,3> vs0(nz,nx,ny);
    for(int i=0;i<ny;i++){
    for(int j=0;j<nx;j++){
    for(int k=0;k<nz;k++){
        vs0(k,j,i) = vs(j,i,k);
    }}}

    // compute dispersion map
    omp_set_num_threads(nthreads);
    #pragma omp parallel for shared(vs0,pv0,mod)
    for(int n=0;n<nx*ny;n++){
        int i = n % nx;
        int j = n/nx;

        // convert vs to vp, rho by empirical relations
        float v[nz],rho[nz],vp[nz];
        for(int k=0;k<nz;k++){
            v[k] = vs0(k,i,j);
            empirical_relation(v+k,vp+k,rho+k);
        }

        // compute frechet kernel
        __dispMap__(&pv0(0,n),v,vp,rho,mod.dep.data(),nz,sublayer,tRc,tRg,tLc,tLg);
    }

    pv = pv0.transpose();
}

/**
 * Compute 1-D surface wave frechet kernel
 * Parameters:
 *  mod     : 3-D initial model
 *  vs      : 3-D shear wave, shape(nx,ny,nz)
 *  pv      : dispersion map, shape(nx*ny,nt)
 *  kernel  : frechet kernel, shape(nx*ny,nz,nt) 
 */
void SurfTime:: SurfWaveKernel(MOD3d &mod,Tensor<float,3> &vs,
                                MatrixXd &pv,Tensor<double,3> &kernel)
{
    // get vs-dimension
    int nx = mod.nx,ny = mod.ny,nz = mod.nz;
    int nt = pv.cols();

    // copy vs to vs0 
    Tensor<float,3> vs0(nz,nx,ny);
    for(int i=0;i<ny;i++){
    for(int j=0;j<nx;j++){
    for(int k=0;k<nz;k++){
        vs0(k,j,i) = vs(j,i,k);
    }}}

    // compute frechet kernel and dispersion map
    MatrixXd pv0(nt,nx*ny);
    Tensor<double,3> kernel0(nt,nz,nx*ny);
    omp_set_num_threads(nthreads);
    #pragma omp parallel for shared(vs0,pv0,mod,kernel0)
    for(int n=0;n<nx*ny;n++){
        int i = n % nx;
        int j = n/nx;

        // convert vs to vp, rho by empirical relations
        float v[nz],rho[nz],vp[nz];
        for(int k=0;k<nz;k++){
            v[k] = vs0(k,i,j);
            empirical_relation(v+k,vp+k,rho+k);
        }
        // compute dispersion map
        __dispMap__(&pv0(0,n),v,vp,rho,mod.dep.data(),nz,sublayer,tRc,tRg,tLc,tLg);

        // compute dispersion kernel
        float dv = 0.01,v0,vp0,rho0;
        double vmap2[nt],vmap1[nt];
        for(int k=0;k<nz;k++){
            v0 = v[k];
            vp0 = vp[k];
            rho0 = rho[k];
            v[k] = v0 * (1.0 + 0.5 * dv);
            empirical_relation(v+k,vp+k,rho+k);
            __dispMap__(vmap2,v,vp,rho,mod.dep.data(),nz,sublayer,tRc,tRg,tLc,tLg);

            v[k] = v0* (1.0 - 0.5 * dv);
            empirical_relation(v+k,vp+k,rho+k);
            __dispMap__(vmap1,v,vp,rho,mod.dep.data(),nz,sublayer,tRc,tRg,tLc,tLg);

            for(int ii=0;ii<nt;ii++) {
                kernel0(ii,k,n) = (vmap2[ii] - vmap1[ii]) / (v0 * dv);
            }
            v[k] = v0;
            vp[k] = vp0;
            rho[k] = rho0;
        }
    }

    // copy pv0 and kernel0 to pv and kernel
    pv = pv0.transpose();
    for(int i=0;i<nt;i++){
    for(int j=0;j<nz;j++){
    for(int k=0;k<nx*ny;k++){
        kernel(k,j,i) = kernel0(i,j,k);
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
    int np = Pairs.size();;
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
                            std::string save_dir)
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
                int nar = (frechet.col(i).array().abs()>0).cast<int>().sum();
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
read_Frechet_Kernel(std::string basedir,csr_matrix<float> &smat)
{
    omp_set_num_threads(nthreads);
    #pragma parallel for shared(smat)
    for(int p=0;p<nthreads;p++){
        std::string filename = basedir + "/" +  std::to_string(p) + ".txt";

        // open file
        char line[100];
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
                int ierr =fscanf(fp,"%d%f\n",smat.indices + i,smat.data + i);
            }
        }

        //close and delete file 
        fclose(fp);
        int ierr = system(("rm -r "+ filename).c_str() );
    } 
}

void SurfTime:: 
write_disper(Eigen::VectorXf &dsyn,std::string filename){
    std::ofstream outfile;
    outfile.open(filename);

    // output
    int c = 0; // counter
    double rad2deg = 180. / pi;
    for(int n=0;n<Pairs.size();n++){
        StationPair &p = Pairs[n];
        float elat = (pi * 0.5 - p.srcx) * rad2deg, elon = p.srcz * rad2deg;
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
            float  slat = (pi*0.5 - p.rcx[i]) * rad2deg, slon = p.rcz[i] * rad2deg;
            outfile << slat << " " << slon << " " << T << " " << p.wavetype << "\n";
            c += 1;
        }
    }
    outfile.close();
}