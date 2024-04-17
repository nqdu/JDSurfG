#include "JSurfGTomo/JSurfGTomo.hpp"
#include "shared/IOFunc.hpp"
#include "shared/gaussian.hpp"
#include "gravity/gravity_module.hpp"
#include "SWD/empirical.hpp"
#include "shared/smoothing.hpp"

void JSurfGTomo::
read_invparams(const std::string &paramfile)
{
    param.read_file(paramfile);

    // read noise level 
    std::ifstream infile; infile.open(paramfile);
    read_par_regex("NOISE_LEVEL_SWD",param.noilevel1,infile);
    read_par_regex("NOISE_LEVEL_GRAV",param.noilevel2,infile);
    read_par_regex("WEIGHT_SWD",param.weight1,infile);
    read_par_regex("WEIGHT_GRAV",param.weight2,infile);
    read_par_regex("RELATIVE_P",param.p,infile);

    if(param.ifsyn == 1) {
        printf("noise level for two dataset: %g %g\n",param.noilevel1,param.noilevel2);
        
    }
    printf("weight factor for two dataset: %g %g\n",param.weight1,param.weight2);
    read_par_regex("REMOVE_AVG",param.rm_grav_avg,infile);

    infile.close();
}

void JSurfGTomo ::
read_model(const std::string &modfile,const std::string &modtrue,
            const std::string &modref)
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
    if(modtrue != "None") {
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

    // read ref model if required
    vsref.resize(nx,ny,nz);
    if(modref != "None") {
        infile.open(modref);
        if(!infile.is_open()) {
            printf("cannot open %s\n",modref.c_str());
            exit(1);
        }
        getline(infile,line); getline(infile,line); getline(infile,line);
        getline(infile,line);
        for(int k=0;k<nz;k++){
        for(int j=0;j<ny;j++){
        for(int i=0;i<nx;i++){
            infile >> vsref(i,j,k);
        }}}
        infile.close();
    } 
    else {
        for(int k = 0; k < nz ;k ++) {
            float mean = 0.;
            for(int j = 0; j < ny -2;j ++){
            for(int i = 0; i < nx -2; i ++){
                mean += vsinit(i+1,j+1,k);
            }}
            mean /= (nx-2) * (ny-2);
            for(int j = 0; j < ny;j ++){
            for(int i = 0; i < nx; i ++){
                vsref(i,j,k) = mean;
            }}
        }
    }
}

void JSurfGTomo:: 
read_data(const std::string &swdfile,const std::string &gravfile)
{
    // read swd data
    surf.read_swd_data(swdfile);

    // read gravity data
    OBSSphGraRandom obssph;
    obssph.read_obs_data(gravfile);

    int np = obssph.np;
    obsg.resize(np); lon_grav.resize(np); lat_grav.resize(np);
    for(int i = 0; i < np; i ++) {
        obsg[i] = obssph.Gr[i];
        lon_grav[i] = obssph.lon[i];
        lat_grav[i] = obssph.lat[i];
    }

    if(param.rm_grav_avg) {
        float mean = obsg.sum() / obsg.size();
        obsg -= mean;
    }
}

void JSurfGTomo:: read_gravmat(const std::string &gravmat)
{
    gmat.read_binary(gravmat);
}

void JSurfGTomo::
compute_gravity(const fmat3 &vs,fvec &dgsyn) const
{
    int m = obsg.size();
    dgsyn.resize(m); dgsyn.setZero();

    // compute gravity
    int nx = vs.dimension(0), ny = vs.dimension(1);
    int nz = vs.dimension(2);
    int size = (nx-2) * (ny - 2) * (nz-1);
    fvec drho(size);
    for(int k = 0; k < nz - 1; k ++) {
    for(int j = 0; j < ny - 2; j ++) {
    for(int i = 0; i < nx - 2; i ++) {
        float b = vs(i+1,j+1,k);
        float rho,rho1,a;
        int n = k * (ny-2) * (nx-2) + j * (nx-2) + i;
        empirical_relation(b,a,rho);
        drho[n] = rho;
        b = vsref(i+1,j+1,k);
        empirical_relation(b,a,rho1);
        drho[n] -= rho1;
    }}}

    // aprod
    gmat.aprod(1,drho.data(),dgsyn.data());

    //remove average value if required
    if(param.rm_grav_avg) {
        dgsyn -= dgsyn.sum() / m;
    }
}

void JSurfGTomo::  
checkerboard()
{
    int m = surf.obst.size();
    fvec dsyn(m);
    surf.travel_time(vstrue,dsyn);

    // add noise
    for(int i=0;i<m;i++){
        //dsyn(i) *= (1.0 + param.noiselevel * gaussian());
        dsyn(i) += (float) (param.noilevel1 * gaussian());
    }
    this -> surf.obst = dsyn*1.0f;

    // gravity
    fvec dgsyn = obsg * 0.;
    this -> compute_gravity(vstrue,dgsyn);
    m = obsg.size();
    for(int i=0;i<m;i++){
        //dsyn(i) *= (1.0 + param.noiselevel * gaussian());
        dgsyn(i) +=  (float) (param.noilevel2 * gaussian());
    }
    this-> obsg = dgsyn * 1.0f;
}

float JSurfGTomo:: 
compute_misfit(const fvec &dsyn) const
{
    // compute weight factor
    int m1 = surf.obst.size(), m2 = this->obsg.size();
    float s1,s2; 
    s1 = param.p / (m1 * std::pow(param.weight1,2));
    s2 = (1. - param.p) / (m1 * std::pow(param.weight2,2));
    s2 = s2 / s1;
    s1 = 1.;

    // misfit
    float chi1 = (dsyn.segment(0,m1) - surf.obst).square().sum() * 0.5;
    float chi2 = (dsyn.segment(m1,m2) - obsg).square().sum() * 0.5;

    return chi1 * s1 + chi2 * s2;
}

void JSurfGTomo:: 
forward(const fvec &x,fvec &dsyn) const
{
    // copy x to vs 
    fmat3 vsf = vsinit;
    int nx = vsf.dimension(0), ny = vsf.dimension(1);
    int nz = vsf.dimension(2);

    int ic = 0;
    for(int k=0;k<nz-1;k++){
    for(int j=0;j<ny-2;j++){
    for(int i=0;i<nx-2;i++){
        vsf(i+1,j+1,k) = x[ic];
        ic += 1;
    }}}

    // allocate space
    int m1 = surf.obst.size(), m2 = this->obsg.size();
    dsyn.resize(m1+m2);

    // forward
    fvec d1,d2;
    surf.travel_time(vsf,d1);
    this -> compute_gravity(vsf,d2);
    dsyn.segment(0,m1) = d1;
    dsyn.segment(m1,m2) = d2;
}

void JSurfGTomo::
compute_grav_grad(const fmat3 &vs,fvec &dgsyn,fvec &grad) const
{
    int m = obsg.size();
    dgsyn.resize(m);

    // compute gravity
    int nx = vs.dimension(0), ny = vs.dimension(1);
    int nz = vs.dimension(2);
    int size = (nx-2) * (ny - 2) * (nz-1);
    fvec drho(size), deriv(size);
    for(int k = 0; k < nz - 1; k ++) {
    for(int j = 0; j < ny - 2; j ++) {
    for(int i = 0; i < nx - 2; i ++) {
        float b = vs(i+1,j+1,k);
        float rho,a;
        float tmp1,tmp2;
        int n = k * (ny-2) * (nx-2) + j * (nx-2) + i;
        empirical_relation(b,a,rho);
        empirical_deriv(a,b,tmp1,tmp2);
        deriv[n] = tmp1 * tmp2;
        drho[n] = rho;
        b = vsref(i+1,j+1,k);
        empirical_relation(b,a,rho);
        drho[n] -= rho;
    }}}

    // synthetic data
    gmat.aprod(1,drho.data(),dgsyn.data()); // dgsyn = G * drho 

    // compute grad
    grad.setZero(); 
    fvec obs = -obsg;
    gmat.aprod(2,grad.data(),dgsyn.data()); // temp = G.T G drho 
    gmat.aprod(2,grad.data(),obs.data()); // temp = -G.T d 

    // multiply vs to grad
    for(int k = 0; k < nz - 1; k ++) {
    for(int j = 0; j < ny - 2; j ++) {
    for(int i = 0; i < nx - 2; i ++) {
        int n = k * (ny-2) * (nx-2) + j * (nx-2) + i;
        grad[n] *= vs(i+1,j+1,k) * deriv[n];
    }}}
}

void JSurfGTomo:: 
compute_grad(const fvec &x,fvec &dsyn,fvec &grad) const
{
    // copy x to vs 
    fmat3 vsf = vsinit;
    int nx = vsf.dimension(0), ny = vsf.dimension(1);
    int nz = vsf.dimension(2);

    int ic = 0;
    for(int k=0;k<nz-1;k++){
    for(int j=0;j<ny-2;j++){
    for(int i=0;i<nx-2;i++){
        vsf(i+1,j+1,k) = x[ic];
        ic += 1;
    }}}

    // allocate space
    int m1 = surf.obst.size(), m2 = this->obsg.size();
    int n = x.size();
    fvec grad1(n),grad2(n);
    fvec d1(m1),d2(m2);
    dsyn.resize(m1+m2);

    // compute swd data/grad
    surf.compute_grad(vsf,d1,grad1);
    this -> compute_grav_grad(vsf,d2,grad2);

    // compute gravity data/grad
    grad2.setZero();
    fvec x1 = x;
    gmat.aprod(1,x1.data(),grad2.data()); // Gm

    // cache data    
    dsyn.segment(0,m1) = d1;
    dsyn.segment(m1,m2) = d2;

    // cache grad
    float s1,s2; 
    s1 = param.p / (m1 * std::pow(param.weight1,2));
    s2 = (1. - param.p) / (m1 * std::pow(param.weight2,2));
    s2 = s2 / s1;
    s1 = 1.;
    grad = grad1 * s1 + grad2 * s2;
}

/**
 * @brief smooth gradient by using gaussian smoothing algorithm
 * 
 * @param grad  gradient array, shape(nlat-2,nlon-2,nz-1), col major
 */
void JSurfGTomo::
smoothing(fvec &grad) const 
{
    int nx = vsinit.dimension(0), ny = vsinit.dimension(1);
    int nz = vsinit.dimension(2);

    // find minimum dz
    float dz = 1000 * (dep[nz-2] - dep[0]); 
    for(int i = 0; i < nz-2; i ++) {
        dz = std::min(dz,dep[i+1] - dep[i]);
    }

    // create regular data
    int nz1 = (dep[nz-2] - dep[0]) / (dz * 0.99);
    dz = (dep[nz-2] - dep[0]) / (nz1 - 1);
    fmat3 gradr(nx-2,ny-2,nz1);
    interp_irregular_z(grad.data(),gradr.data(),nx-2,ny-2,nz-1,nz1,dep.data(),true);

    // smoothing
    if (param.smooth_in_km) {
        float dx = std::abs(lat[1] - lat[0]);
        float dy = std::abs(lon[1] - lon[0]);
        smooth_sph_pde(gradr.data(),nx-2,ny-2,nz1,dx,dy,dz,
                        lat[1],lon[1],dep[0],param.sigma_h,
                        param.sigma_v);
    }
    else {
        smooth_cart_pde(gradr.data(),nx-2,ny-2,nz1,param.sigma_h,param.sigma_v);
    }

    // interpolate back
    interp_irregular_z(grad.data(),gradr.data(),nx-2,ny-2,nz-1,nz1,dep.data(),false);
}

/**
 * @brief assemble gravity and SWD matrix to a global one 
 * @param vsf current model
 * @param smat global frechet matrix
 * @param gmat gravity frechet matrix
 * @param m1/m2 # of data for swd/gravity
 * @param s1/s2 weight factor for swd/gravity
*/
static void 
assemble(const fmat3 &vsf,const csr_matrix &gmat,csr_matrix &smat,
        int m1,int m2,float s1,float s2)
{
    // get dimension
    int nx = vsf.dimension(0), ny = vsf.dimension(1);

    // loop around all rows to update index
    for(int r = 0; r < m2 ;r ++){
        int start = gmat.indptr[r];
        int end = gmat.indptr[r + 1];
        int rwc = m1 + r;
        smat.indptr[rwc + 1]  = smat.indptr[rwc] + end - start;
    }

    // assemble
    #pragma omp parallel for shared(smat,vsf,gmat)
    for(int r = 0; r < m2; r ++){
        int start = gmat.indptr[r];
        int end = gmat.indptr[r + 1];
        int rwc = m1 + r;

        for(int c = start; c < end; c ++){
            int col = gmat.indices[c];
            int k = col / ((nx-2) * (ny-2));
            int ridx = col % ((nx-2) * (ny-2));
            int i = ridx %(nx-2), j = ridx / (nx-2);

            // compute derivaties
            float vp,vs,rho;
            vs = vsf(i+1,j+1,k);
            empirical_relation(vs,vp,rho);
            float tmp1,tmp2;
            empirical_deriv(vp,vs,tmp1,tmp2);

            // append gmat to smat
            int count = smat.indptr[rwc] + c - start;
            smat.data[count] = gmat.data[c] * tmp1 * tmp2;
            smat.indices[count] = gmat.indices[c];
        }
    }

    // add weights
    for(int i = 0; i < m1 + m2; i ++){
    for(int j = smat.indptr[i]; j < smat.indptr[i+1]; j ++){
        float w =  i < m1 ? s1 : s2;
        smat.data[j] *= w;
    }}
}

void JSurfGTomo::
inversion(fmat3 &vsf,fvec &dsyn) const
{
    int m1 = surf.obst.size(), m2 = obsg.size();
    int nx = vsf.dimension(0), ny = vsf.dimension(1);
    int nz = vsf.dimension(2);  
    int n = (nx-2) * (ny -2) * (nz  - 1);

    // residual vector and velocity variation vector
    fvec res(m1 + m2 + n),dv(n);
    dsyn.resize(m1+m2); 
    fvec dsyn1(m1),dsyn2(m2);
    res.setZero();
    dv.setZero();

    // comptue gravity
    this -> compute_gravity(vsf,dsyn2);

    // compute and save frechet matrix
    int nar = surf.frechet_matrix(vsf,dsyn1,"frechet.bin");
    res.segment(0,m1) = surf.obst - dsyn1;
    res.segment(m1,m2) = obsg - dsyn2;

    // save dsyn
    dsyn.segment(0,m1) = dsyn1;
    dsyn.segment(m1,m2) = dsyn2;

    // compute weights
    float s1,s2;
    if(param.ifsyn) {
        s1 = param.noilevel1 / std::sqrt(m1 * 1.);
        s2 = param.noilevel2 / std::sqrt(m1 * 1.);
    }
    else {
        s1 = param.weight1 / std::sqrt(m1 * 1.);
        s2 = param.weight2 / std::sqrt(m1 * 1.);
    }
    s1 = std::sqrt(param.p / m1) / s1;
    s2 = std::sqrt((1-param.p) / m2) / s2;
    //s2 = s2 / s1; s1 = 1;

    // initialize matrix and read frechet binary file
    csr_matrix smat(m1+m2+n,n,nar + gmat.nonzeros + n*7);
    FILE *fp = fopen("frechet.bin","rb");
    assert(fread(smat.indptr,sizeof(int),1,fp) == 1);
    assert(fread(smat.indptr,sizeof(int),1,fp) == 1);
    smat.indptr[0] = 0;
    for(int i=0;i<m1;i++){
        int col,nar1;
        assert(fread(&col,sizeof(int),1,fp) == 1);
        assert(fread(&nar1,sizeof(int),1,fp) == 1);
        int start = smat.indptr[col];
        smat.indptr[col+1] = start +  nar1;

        // read data
        assert(fread(smat.indices+start,sizeof(int),nar1,fp) == (size_t)nar1);
        assert(fread(smat.data+start,sizeof(float),nar1,fp) == (size_t)nar1);
    }
    fclose(fp);
    for(int i=m1;i<m1+m2+n;i++) smat.indptr[i+1] = smat.indptr[m1];

    // delete kernel file
    std::remove("frechet.bin");

    // add weights to res 
    for(int i = 0; i < m1; i ++) res[i] *= s1;
    for(int i = m1; i < m1+m2; i ++) res[i] *= s2;

    // assemble
    assemble(vsf,gmat,smat,m1,m2,s1,s2);

    // add regularization
    add_regularization(smat,param.smooth,nx-2,ny-2,nz-1);

    // solve equations by lsmr
    printf("solving linear systems by LSMR ...\n");
    int itnlim = n * 2; 
    LSMRDict<float> dict(itnlim,10,param.damp,param.smooth,param.nthreads);
    dict.verbose = false; 
    smat.lsmr_solver(dv.data(),res.data(),dict);
    printf("max negative and positive perturbation:  %f %f\n",
            dv.minCoeff(),dv.maxCoeff());
    
    // tackle large variations
    Eigen::TensorMap<fmat3>dx(dv.data(),nx-2,ny-2,nz-1);
    float minvel = param.minvel;
    float maxvel = param.maxvel;
    for(int k=0;k<nz-1;k++){
    for(int j=0;j<ny-2;j++){
    for(int i=0;i<nx-2;i++){
        // if(dx(i,j,k) > 0.5) dx(i,j,k) = 0.5;
        // if(dx(i,j,k) < -0.5) dx(i,j,k) = -0.5;
        float temp = dx(i,j,k) + vsf(i+1,j+1,k);
        if(temp > maxvel) temp = maxvel;
        if(temp < minvel) temp = minvel;
        if(vsf(i+1,j+1,k) !=0.0) vsf(i+1,j+1,k) = temp;
    }}}
}