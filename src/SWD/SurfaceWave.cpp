#include "SWD/SurfaceWave.hpp"
#include "SWD/swd.hpp"
#include "shared/parallel.hpp"
#include "SWD/empirical.hpp"
#include "shared/csr_matrix.hpp"
#include<fstream>

/**
 * compute global index in [0,kmax-1] for a given wave type and given mode
 * @param period_idx period index for this wave type
 * @param mode mode =0 fundamental, = 1 first order
 * @param wavetype wavetype, one of [Rc,Rg,Lc,Lg]
 */ 
int SurfTime :: 
get_period_index(int period_idx,int mode,const std::string &wavetype) const
{
    int idx;
    if(wavetype == "Rc") {
        int kmax = 0;
        for(int i=0; i < mode;i++){
            kmax += kmaxRc[i];
        }
        idx = period_idx + kmax;
    }
    else if (wavetype == "Rg") {
        int kmax = kmaxRc.sum();
        for(int i=0; i < mode;i++){
            kmax += kmaxRg[i];
        }
        idx = period_idx + kmax;
    }
    else if (wavetype == "Lc") {
        int kmax = kmaxRc.sum() + kmaxRg.sum();
        for(int i=0; i < mode;i++){
            kmax += kmaxLc[i];
        }
        idx = period_idx + kmax; 
    }
    else {
        int kmax = kmaxRc.sum() + kmaxRg.sum() + kmaxLc.sum();
        for(int i=0; i < mode;i++){
            kmax += kmaxLg[i];
        }
        idx = period_idx + kmax; 
    }

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
static void 
_ParamConvert(const float *dep,const float *vp,const float *vs,
             const float *rho,int nz,int sublayer,float *rthk,
             float *rvp,float *rvs,float *rrho,int rmax,
             double *rdcdm,double *dcdm)
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
 * refine grid based model to layerd based model, water layers considered
 * 
 * \param dep,vp,vs,rho grid-based model
 * \param nz no. of points in grid-based model
 * \param sublayer insert sublayer points to construct layered model
 * \param rdep, rvp, rvs, rrho, rthk layered velocity model 
 * 
 * \return rmax no. of layers of refined model
 */
static int 
grid2LayerModel(const float *dep,const float *vp,const float *vs,
                const float *rho,int nz,int sublayers,
                float *rthk,float *rvp,float *rvs,float *rrho)
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

static int get_nt_row(const dvec &tRc) {
    int n;
    int col = tRc.size();
    for(n = 0; n < col; n ++){
        if(tRc[n] == -1){
            break;
        }
    }

    return n;
}

/**
 * internal function to compute disperion for one cell
 * @param  vs,vp,rho,dep parameters for one cell
 * @param  nz  no. of points in vertical direction
 * @param sublayer insert sublayer layers in adjacent points
 * @param tRc... period vector for each wavetype/mode
 * @param vc,vout dispersion for one cell, size(kmaxRc+kmaxRg+kmaxLc+kmaxLg)
*/
static void
disp1D(const float *vs,const float *vp,const float *rho,const float *dep,
        int nz,const dmat2 &tRc,const dmat2 &tRg,
        const dmat2 &tLc,const dmat2 &tLg,double *vc,double *vout)
{
    // convert grid-based model to layer-based model
    const int sublayer = 3;
    int rmax = nz + (nz-1) * sublayer;
    float *rvp = new float [rmax],*rvs = new float [rmax];
    float *rthk = new float [rmax],*rrho = new float [rmax];
    rmax = grid2LayerModel(dep,vp,vs,rho,nz,sublayer,rthk,rvp,rvs,rrho);
    
    // prepare parameters
    int kmaxRc = tRc.size(),kmaxRg = tRg.size();
    int kmaxLc= tLc.size(), kmaxLg = tLg.size();
    bool sphere = true, keep_flat = false;

    // lambda function to find the 

    // compute dispersion
    int idx = 0;
    if(kmaxRc > 0){
        int nmode = tRc.cols();
        for(int mode = 0; mode < nmode; mode ++) {
            int nt = get_nt_row(tRc.col(mode));
            if(nt > 0) {
                surfdisp(rthk,rvp,rvs,rrho,rmax,&tRc(0,mode),vc+idx,
                        nt,"Rc",mode,sphere,keep_flat);
                memcpy(vout+idx,vc+idx,nt*sizeof(double));
                idx += nt;
            }
        }
    }
    if(kmaxRg > 0){
        int nmode = tRg.cols();
        for(int mode = 0; mode < nmode; mode ++) {
            int nt = get_nt_row(tRg.col(mode));
            if(nt > 0) {
                surfdisp(rthk,rvp,rvs,rrho,rmax,&tRg(0,mode),vc+idx,
                        nt,"Rc",mode,sphere,keep_flat);
                groupvel_r(rthk,rvp,rvs,rrho,rmax,&tRg(0,mode),
                            vout+ idx,nt,mode,sphere);
                idx += nt;
            }
        }
    }
    if(kmaxLc > 0){
        int nmode = tLc.cols();
        for(int mode = 0; mode < nmode; mode ++) {
            int nt = get_nt_row(tLc.col(mode));
            if(nt > 0) {
                surfdisp(rthk,rvp,rvs,rrho,rmax,&tLc(0,mode),vc+idx,
                        nt,"Lc",mode,sphere,keep_flat);
                memcpy(vout+idx,vc+idx,nt*sizeof(double));
                idx += nt;
            }
        }
    }

    if(kmaxLg > 0){
        int nmode = tLg.cols();
        for(int mode = 0; mode < nmode; mode ++) {
            int nt = get_nt_row(tLg.col(mode));
            if(nt > 0) {
                surfdisp(rthk,rvp,rvs,rrho,rmax,&tLg(0,mode),vc+idx,
                        nt,"Lc",mode,sphere,keep_flat);
                groupvel_l(rthk,rvs,rrho,rmax,&tLg(0,mode),
                            vout + idx,nt,mode,sphere);
                idx += nt;
            }
        }
    }

    // free space
    delete[] rrho;
    delete[] rvp;
    delete[] rvs;
    delete[] rthk;
}

/**
 * internal function to compute disperion for one cell
 * @param  vs,vp,rho,dep parameters for one cell
 * @param  nz  no. of points in vertical direction
 * @param sublayer insert sublayer layers in adjacent points
 * @param tRc... period vector for each wavetype/mode
 * @param vc,vout dispersion for one cell, size(kmaxRc+kmaxRg+kmaxLc+kmaxLg)
 * @param kvs/vp/rho sensitivity kernel for vs/vp/rho, shape(nz,nt),column major
*/
static void 
kernel1D(const float *vs,const float *vp,const float *rho,const float *dep,int nz,
         const dmat2 &tRc,const dmat2 &tRg,const dmat2 &tLc,
         const dmat2 &tLg,double *vc,double *vout,double *kvs,double *kvp,
         double *krho)
{
    // convert grid-based model to layer-based model
    const int sublayer = 3;
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
    int iflsph = 1;

    // compute dispersion and sensitivity kernel
    int idx = 0;
    if(kmaxRc > 0){ // rayleigh phase
        int nmode = tRc.cols();
        for(int mode = 0; mode < nmode; mode ++) {
            int nt = get_nt_row(tRc.col(mode));
            if(nt > 0) {
                surfdisp(rthk,rvp,rvs,rrho,rmax,&tRc(0,mode),vc+idx,
                        nt,"Rc",mode,sphere);
                for(int i = 0; i < nt; i ++) {
                    int k = i + idx;
                    double cg;
                    sregn96_(rthk,rvp,rvs,rrho,rmax,&tRc(i,mode),vc + k,&cg,
                            ur,uz,tr,tz,dcdar,dcdbr,dcdhr,dcdrr,iflsph);
                    _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                                rthk,rvp,rvs,rrho,rmax,dcdar,kvp+k*nz);
                    _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                                rthk,rvp,rvs,rrho,rmax,dcdbr,kvs+k*nz);
                    _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                                rthk,rvp,rvs,rrho,rmax,dcdrr,krho+k*nz);
                    vout[k] = vc[k];
                }
                idx += nt;
            }
        } 
    }
    if(kmaxRg > 0){
        int nmode = tRg.cols();
        for(int mode = 0; mode < nmode; mode ++) {
            int nt = get_nt_row(tRg.col(mode));
            double cp[nt],t1[nt],t2[nt],c1[nt],c2[nt];
            double dt = 0.01;
            for(int i=0;i<nt;i++){
                t1[i] = tRg(i,mode) * (1.0 + 0.5 * dt);
                t2[i] = tRg(i,mode) * (1.0 - 0.5 * dt);
            }
            surfdisp(rthk,rvp,rvs,rrho,rmax,&tRg(0,mode),cp,nt,"Rc",mode,sphere);
            surfdisp(rthk,rvp,rvs,rrho,rmax,t1,c1,nt,"Rc",mode,sphere);
            surfdisp(rthk,rvp,rvs,rrho,rmax,t2,c2,nt,"Rc",mode,sphere);
            for(int i=0;i<nt;i++){
                int k = i + idx;
                sregnpu_(rthk,rvp,rvs,rrho,rmax,&tRg(i,mode),cp+i,vc+k,
                        ur,uz,tr,tz,t1+i,c1+i,t2+i,c2+i,dcdar,
                        dcdbr,dcdhr,dcdrr,dudar,dudbr,dudhr,dudrr,iflsph);
                _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                            rthk,rvp,rvs,rrho,rmax,dudar,kvp+k*nz);
                _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                            rthk,rvp,rvs,rrho,rmax,dudbr,kvs+k*nz);
                _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                            rthk,rvp,rvs,rrho,rmax,dudrr,krho+k*nz);
                vout[k] = vc[k];
                vc[k] = cp[i];        
            } 
            idx += nt;
        }
    }
    if(kmaxLc > 0){ // Love phase
        int nmode = tLc.cols();
        for(int mode = 0; mode < nmode; mode ++) {
            int nt = get_nt_row(tLc.col(mode));
            if(nt > 0) {
                surfdisp(rthk,rvp,rvs,rrho,rmax,&tLc(0,mode),vc+idx,
                        nt,"Lc",mode,sphere);
                for(int i = 0; i < nt; i ++) {
                    int k = i + idx;
                    double cg;
                    sregn96_(rthk,rvp,rvs,rrho,rmax,&tLc(i,mode),vc + k,&cg,
                            ur,uz,tr,tz,dcdar,dcdbr,dcdhr,dcdrr,iflsph);
                    _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                                rthk,rvp,rvs,rrho,rmax,dcdar,kvp+k*nz);
                    _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                                rthk,rvp,rvs,rrho,rmax,dcdbr,kvs+k*nz);
                    _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                                rthk,rvp,rvs,rrho,rmax,dcdrr,krho+k*nz);  
                    vout[k] = vc[k];
                }
                idx += nt;
            }
        } 
    }
    if(kmaxLg > 0){
        int nmode = tLg.cols();
        for(int mode = 0; mode < nmode; mode ++) {
            int nt = get_nt_row(tLg.col(mode));
            double cp[nt],t1[nt],t2[nt],c1[nt],c2[nt];
            double dt = 0.01;
            for(int i=0;i<nt;i++){
                t1[i] = tLg(i,mode) * (1.0 + 0.5 * dt);
                t2[i] = tLg(i,mode) * (1.0 - 0.5 * dt);
            }
            surfdisp(rthk,rvp,rvs,rrho,rmax,&tLg(0,mode),cp,nt,"Lc",mode,sphere);
            surfdisp(rthk,rvp,rvs,rrho,rmax,t1,c1,nt,"Lc",mode,sphere);
            surfdisp(rthk,rvp,rvs,rrho,rmax,t2,c2,nt,"Lc",mode,sphere);
            for(int i=0;i<nt;i++){
                int k = i + idx;
                sregnpu_(rthk,rvp,rvs,rrho,rmax,&tLg(i,mode),cp+i,vc+k,
                        ur,uz,tr,tz,t1+i,c1+i,t2+i,c2+i,dcdar,
                        dcdbr,dcdhr,dcdrr,dudar,dudbr,dudhr,dudrr,iflsph);
                _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                            rthk,rvp,rvs,rrho,rmax,dudar,kvp+k*nz);
                _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                            rthk,rvp,rvs,rrho,rmax,dudbr,kvs+k*nz);
                _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                            rthk,rvp,rvs,rrho,rmax,dudrr,krho+k*nz);
                vout[k] = vc[k];
                vc[k] = cp[i];     
            } 
            idx += nt;
        }
    }

    // free space
    delete[] rrho;delete[] rvp;delete[] rvs;delete[] rthk;
    delete[] ur; delete [] uz; delete [] tr; delete[] tz;
    delete [] dudrr; delete[] dudar; delete[] dudbr; delete[] dudhr;
    delete [] dcdrr; delete[] dcdar; delete[] dcdbr; delete[] dcdhr;
}

/**
 * Compute 2-Ddispersion map
 * @param vs vs velocity (nx,ny,nz)
 * @param vc phase dispersion map(nx*ny,nt)
 * @param vout used dispersion map(nx*ny,nt),
*/
void SurfTime ::
get_2d_map(const fmat3 &vs,dmat2 &vc,dmat2 &vout) const
{
    // get dimension
    int nx = vs.dimension(0),ny = vs.dimension(1),nz = vs.dimension(2);
    int nt = kmax;

    // copy vs to vs0
    vc.resize(nx*ny,nt); vout.resize(nx*ny,nt);

    // compute dispersion map
    #pragma omp parallel for shared(vs,vc,vout)
    for(int n=0;n<nx*ny;n++){
        int i = n % nx;
        int j = n/nx;

        // convert vs to vp, rho by empirical relations
        float v[nz],rho[nz],vp[nz],dep[nz];
        for(int k=0;k<nz;k++){
            v[k] = vs(i,j,k);
            dep[k] = this->depth[k];
            empirical_relation(v[k],vp[k],rho[k]);
        }

        // compute frechet kernel
        dvec c(nt),cout(nt);
        disp1D(v,vp,rho,dep,nz,tRc,tRg,tLc,tLg,c.data(),cout.data());
        vc.row(n) = c.transpose();
        vout.row(n) = cout.transpose();
    }
}

/**
 * Compute 1-D surface wave frechet kernel
 *  @param vs 3-D shear wave, shape(nx,ny,nz)
 *  @param vmap dispersion map, shape(nx*ny,nt)
 *  @param kernel frechet kernel, shape(nx*ny,nz,nt) 
 */
void SurfTime:: 
get_1d_kernel(const fmat3 &vs,dmat2 &vc,dmat2 &vout,dmat3 &kernel ) const
{
    // get dimension
    int nx = vs.dimension(0),ny = vs.dimension(1),nz = vs.dimension(2);

    int nt = kmax;
    vc.resize(nx*ny,nt); vout.resize(nx*ny,nt);

    // compute frechet kernel and dispersion map
    kernel.resize(nx*ny,nz,nt);
    #pragma omp parallel for shared(vs,vc,vout,kernel)
    for(int n=0;n<nx*ny;n++){
        int i = n % nx;
        int j = n/nx;

        // convert vs to vp, rho by empirical relations
        float v[nz],rho[nz],vp[nz],dep[nz];
        for(int k=0;k<nz;k++){
            v[k] = vs(i,j,k);
            dep[k] = this -> depth[k];
            empirical_relation(v[k],vp[k],rho[k]);
        }

        // compute dispersion map for vs,vp and rho
        dmat2 kvs(nz,nt),kvp(nz,nt),krho(nz,nt);
        dvec cp(nt),cg(nt);
        kvs.setZero(); kvp.setZero();krho.setZero();
        kernel1D(v,vp,rho,dep,nz,tRc,tRg,tLc,tLg,cp.data(),cg.data(),kvs.data(),
                kvp.data(),krho.data());
        vc.row(n) = cp.transpose();
        vout.row(n) = cg.transpose();

        // convert to vs by empirical relations
        for(int ii=0;ii<nt;ii++){
            float drda,dadb;
            for(int jj=0;jj<nz;jj++){
                if(v[jj] == 0.0 && jj == 0) continue; // water layer
                empirical_deriv(vp[jj],v[jj],drda,dadb);
                kernel(n,jj,ii) = kvs(jj,ii) + kvp(jj,ii) * dadb 
                                + krho(jj,ii) * drda * dadb;
            }
        }
    }    
}

/**
 * compute travel time for each station pair
 * @param vs vs model, shape(nx,ny,nz)
 * @param data shape(nt)
*/
void SurfTime:: 
travel_time(const fmat3 &vs,fvec &data) const
{
    // get vs-dimension
    int nx = vs.dimension(0),ny = vs.dimension(1);
    dmat2 vc,vout; // shape(nx*ny,kmax);
    data.resize(obst.size());

    // compute dispersion map
    this -> get_2d_map(vs,vc,vout);
    // get traveltime for every source-receiver pair
    int np = Pairs.size();

    #pragma omp parallel for shared(vs,vc,vout,data)
    for(int n=0;n<np;n++){
        const StationPair &p = Pairs[n];
        int nr = p.nr;
        int idx = get_period_index(p.period_idx,p.mode,p.wavetype);
        float ttime[nr];

        // run fmst
        fmst_init(nx,ny,lat0,lon0,dlat,dlon);
        fmst_run(&vc(0,idx),p.srcx,p.srcz,p.rcx.data(),p.rcz.data(),nr,ttime);
        fmst_reset(&vout(0,idx));

        for(int i=0;i<nr;i++){
            if(p.wavetype == "Rg" || p.wavetype == "Lg") {
                float fdm[nx*ny];
                ttime[i] = fmst_raypath(p.srcx,p.srcz,p.rcx[i],p.rcz[i],fdm);
            }
            data[p.counter + i] = ttime[i];
        }

        // destroy mem
        fmst_finalize();
    }
}

/**
 * compute frechet matrix 
 * @param vs input vs model, shape(nx,ny,nz)
 * @param data synthetic data vector
 * @param outfile string, output matrix file
 * @return nonzeros  nonzeros in frechet matrix
*/
int SurfTime::
frechet_matrix(const fmat3 &vs,fvec &data,const std::string &outfile) const
{
    // get dimension
    int nx = vs.dimension(0),ny = vs.dimension(1),nz = vs.dimension(2);
    dmat2 vc(nx*ny,kmax),vout(nx*ny,kmax);
    dmat3 kernel(nx*ny,nz,kmax);

    //compute disper map and 1-D kernel
    printf("computing Surface Wave 1-D Frechet Kernel ...\n");
    this -> get_1d_kernel(vs,vc,vout,kernel);

    // save a dispermap 
    // FILE *f1 = fopen("phase.txt","w");
    // for(int i = 0; i < nx*ny;i++){
    //     fprintf(f1,"%g\n",vc(i,22));
    // }
    // fclose(f1);
    // f1 = fopen("group.txt","w");
    // for(int i = 0; i < nx*ny;i++){
    //     fprintf(f1,"%g\n",vout(i,22));
    // }
    // fclose(f1);

    // get omp threads
    int nprocs{};
    #pragma omp parallel
    {
        nprocs = omp_get_num_threads();
    }
    ivec nonzeros(nprocs); nonzeros.setZero();

    // compute traveltime/2D kernel
    printf("computing traveltimes and 2-D Frechet Kernel ...\n");
    int np = Pairs.size(); // no. of pairs
    #pragma omp parallel for shared(vs,vc,vout,data,kernel)
    for(int myrank = 0; myrank < nprocs; myrank ++) {
        // get matrix row/col
        int n = (nx-2)*(ny-2)*(nz-1);
        int m = obst.size();

        // open a file to save kernel
        std::string filename = outfile +  "."  +  std::to_string(myrank);
        FILE *fp;
        if((fp=fopen(filename.c_str(),"wb"))==NULL){
            printf("cannot open file  %s\n",filename.c_str());
            exit(1);
        }
        fwrite(&m,sizeof(int),1,fp);
        fwrite(&n,sizeof(int),1,fp);

        // alloc tasks
        int starid,endid;
        allocate_tasks(np,nprocs,myrank,starid,endid);

        // loop around all station pairs
        for(int ievt = starid; ievt <= endid; ievt ++) {
            const StationPair &p = Pairs[ievt];
            int nr = p.nr;
            int idx = get_period_index(p.period_idx,p.mode,p.wavetype);
            float ttime[nr];

            // compute 3-D pseudo sensitivity kernel
            fmst_init(nx,ny,lat0,lon0,dlat,dlon);
            fmst_run(&vc(0,idx),p.srcx,p.srcz,&p.rcx[0],&p.rcz[0],nr,ttime);
            fmst_reset(&vout(0,idx));

            // loop around all receivers
            for(int i=0;i<nr;i++){
                fmat2 fdm(nx,ny);
                float t = fmst_raypath(p.srcx,p.srcz,p.rcx[i],p.rcz[i],fdm.data());
                if(p.wavetype == "Rg" || p.wavetype == "Lg") {
                    ttime[i] = t;
                }

                // compute frechet kernel in 3d
                fvec frechet(n); frechet.setZero();
                int c = 0;
                for(int iz = 0; iz < nz-1; iz ++ ) {
                for(int iy = 0; iy < ny-2; iy ++ ) {
                for(int ix = 0; ix < nx-2; ix ++ ) {
                    int ii = (iy+1) * nx + ix + 1;
                    frechet(c) = fdm(ix+1,iy+1) * kernel(ii,iz,idx); 
                    c += 1;
                }}}

                int nar = (frechet.abs()>0.0).cast<int>().sum();
                c = p.counter + i;
                fwrite(&c,sizeof(int),1,fp);
                fwrite(&nar,sizeof(int),1,fp);
                nonzeros[myrank] += nar;

                // save kernel
                int myindx[nar]; 
                float value[nar];
                int counter = 0;
                for(int j=0;j<n;j++){
                    float tmp = frechet[j];
                    if(std::abs(tmp)> 0.0) {
                        myindx[counter] = j;
                        value[counter] = tmp;
                        counter += 1;
                    }
                }
                fwrite(myindx,sizeof(int),nar,fp);
                fwrite(value,sizeof(float),nar,fp);
                data[p.counter + i] = ttime[i];
            }

            // destroy memory
            fmst_finalize();
        }

        // close file
        fclose(fp);
    }

    // merge all files to 1 big file
    merge_csr_files(nprocs,outfile);

    // return nonzeros
    return nonzeros.sum();
}

void SurfTime :: 
set_model(const fvec &dep,float goxd,float gozd,float dvxd,float dvzd)
{
    this -> depth = dep;
    lat0 = goxd;
    lon0 = gozd;
    dlat = dvxd;
    dlon = dvzd;
}

void SurfTime:: 
compute_grad(const fmat3 &vs,fvec &data,fvec &grad) const
{
    // get dimension
    int nx = vs.dimension(0),ny = vs.dimension(1),nz = vs.dimension(2);
    dmat2 vc(nx*ny,kmax),vout(nx*ny,kmax);
    dmat3 kernel(nx*ny,nz,kmax);

    //compute disper map and 1-D kernel
    printf("computing SWD time gradient ...\n");
    this -> get_1d_kernel(vs,vc,vout,kernel);

    // get omp threads
    int nprocs{};
    #pragma omp parallel
    {
        nprocs = omp_get_num_threads();
    }

    // compute traveltime/2D kernel
    //printf("computing traveltimes and 2-D Frechet Kernel ...\n");

    // set grad to zero
    int n = (nx-2)*(ny-2)*(nz-1);
    int m = obst.size();
    int np = Pairs.size();
    data.resize(m);
    fmat2 grad_all; grad_all.resize(n,nprocs);
    grad_all.setZero();
 
    // parallelization
    #pragma omp parallel for shared(vs,vc,vout,data,kernel,grad_all)
    for(int myrank = 0; myrank < nprocs; myrank ++) {
        // alloc tasks
        int starid,endid;
        allocate_tasks(np,nprocs,myrank,starid,endid);

        // loop around all station pairs
        for(int ievt = starid; ievt <=endid; ievt ++){
            const StationPair &p = Pairs[ievt];
            int nr = p.nr;
            int idx = get_period_index(p.period_idx,p.mode,p.wavetype);
            float ttime[nr];

            // compute 3-D pseudo sensitivity kernel
            fmst_init(nx,ny,lat0,lon0,dlat,dlon);
            fmst_run(&vc(0,idx),p.srcx,p.srcz,&p.rcx[0],&p.rcz[0],nr,ttime);
            fmst_reset(&vout(0,idx));
            
            // loop all receivers to add grads
            for(int ir = 0; ir < nr; ir ++) {
                fmat2 fdm(nx,ny);
                float t = fmst_raypath(p.srcx,p.srcz,p.rcx[ir],p.rcz[ir],fdm.data());
                if(p.wavetype == "Rg" || p.wavetype == "Lg") {
                    ttime[ir] = t;
                }

                // compute frechet kernel in 3d
                fvec frechet(n); frechet.setZero();
                int c = 0;
                for(int iz = 0; iz < nz-1; iz ++ ) {
                for(int iy = 0; iy < ny-2; iy ++ ) {
                for(int ix = 0; ix < nx-2; ix ++ ) {
                    int ii = (iy+1) * nx + ix + 1;
                    frechet(c) = fdm(ix+1,iy+1) * kernel(ii,iz,idx);  
                    c += 1;
                }}}
                data(p.counter + ir) = ttime[ir];
                grad_all.col(myrank) += frechet * (ttime[ir] - obst(p.counter + ir));
            }

            // destroy mem
            fmst_finalize();
        }
    }

    // add contributations to grad
    grad = grad_all.rowwise().sum();

    // we invert the relative change, so we times vs to grad
    for(int k=0;k<nz-1;k++){
    for(int j=0;j<ny-2;j++){
    for(int i=0;i<nx-2;i++){
        int ic = k * (ny-2) * (nx-2) + j * (nx-2) + i;
        grad(ic) *= vs(i+1,j+1,k);
    }}}
}


void SurfTime:: 
write_syn(const fvec &dsyn, const std::string &outfile) const{

    // open file
    FILE *fp = fopen(outfile.c_str(),"w");

    // output
    int c = 0; // counter
    const double rad2deg = 180. / M_PI;
    for(const auto &p : Pairs){
        float elat = (M_PI * 0.5 - p.srcx) * rad2deg, elon = p.srcz * rad2deg;
        float T;
        if(p.wavetype == "Rc")
            T = tRc(p.period_idx,p.mode);
        else if(p.wavetype == "Rg")
            T = tRg(p.period_idx,p.mode);
        else if(p.wavetype == "Lc")
            T = tLc(p.period_idx,p.mode);
        else 
            T = tLg(p.period_idx,p.mode);
        for(int i=0;i<p.nr;i++){
            float  slat = (M_PI*0.5 - p.rcx[i]) * rad2deg, slon = p.rcz[i] * rad2deg;
            fprintf(fp,"%f %f %f ",sta_dist[c],obst[c],dsyn[c]);
            fprintf(fp,"%f %f %f %f ",elat,elon,slat,slon);
            fprintf(fp,"%f %d %s\n",T,p.mode,p.wavetype.c_str());
            c += 1;
        }
    }

    // close file
    fclose(fp);
}