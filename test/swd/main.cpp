#include "SWD/swd.hpp"
#include <iostream>
#include "numerical.hpp"

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

void compute_kernel(const float *vs,const float *vp,const float *rho,const float *dep,int nz,
                    const double *t,int nt,double *cout,double *kvs,double *kvp,
                    double *krho,std::string &wavetype)
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
    bool sphere = true;
    int iflsph = 1;
    int mode = 0;

    // compute dispersion and sensitivity kernel
    int idx = 0;
    if(wavetype == "Rc"){ // rayleigh phase
        surfdisp(rthk,rvp,rvs,rrho,rmax,t,cout,nt,"Rc",mode,sphere);
        for(int i = 0; i < nt ;i ++) {
            double cg;
            int k = i;
            sregn96_(rthk,rvp,rvs,rrho,rmax,t + i,cout+i,&cg,
                    ur,uz,tr,tz,dcdar,dcdbr,dcdhr,dcdrr,iflsph);
            _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dcdar,kvp+k*nz);
            _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dcdbr,kvs+k*nz);
            _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dcdrr,krho+k*nz);
        }
    }

    if(wavetype == "Rg") {
        double cp[nt],t1[nt],t2[nt],c1[nt],c2[nt];
        double dt = 0.01;
        for(int i=0;i<nt;i++){
            t1[i] = t[i]* (1.0 + 0.5 * dt);
            t2[i] = t[i]* (1.0 - 0.5 * dt);
        }
        surfdisp(rthk,rvp,rvs,rrho,rmax,t,cp,nt,"Rc",mode,sphere);
        surfdisp(rthk,rvp,rvs,rrho,rmax,t1,c1,nt,"Rc",mode,sphere);
        surfdisp(rthk,rvp,rvs,rrho,rmax,t2,c2,nt,"Rc",mode,sphere);
        for(int i=0;i<nt;i++){
            int k = i;
            sregnpu_(rthk,rvp,rvs,rrho,rmax,t + i,cp+i,cout+k,
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


    if(wavetype == "Lc"){ // Love phase
        surfdisp(rthk,rvp,rvs,rrho,rmax,t,cout,nt,"Lc",mode,sphere);
        for(int i = 0; i < nt; i ++) {
            int k = i;
            double cg;
            slegn96_(rthk,rvs,rrho,rmax,t + i,cout+k,&cg,
                    uz,tz,dcdbr,dcdhr,dcdrr,iflsph);
            _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dcdbr,kvs+k*nz);
            _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dcdrr,krho+k*nz);

            // zero out kvp
            for(int iz = 0; iz < nz; iz ++) kvp[k * nz + iz] = 0.;
        }
    }

    if(wavetype == "Lg"){
        double cp[nt],t1[nt],t2[nt],c1[nt],c2[nt];
        double dt = 0.01;
        for(int i=0;i<nt;i++){
            t1[i] = t[i] * (1.0 + 0.5 * dt);
            t2[i] = t[i] * (1.0 - 0.5 * dt);
        }
        surfdisp(rthk,rvp,rvs,rrho,rmax,t,cp,nt,"Lc",mode,sphere);
        surfdisp(rthk,rvp,rvs,rrho,rmax,t1,c1,nt,"Lc",mode,sphere);
        surfdisp(rthk,rvp,rvs,rrho,rmax,t2,c2,nt,"Lc",mode,sphere);
        for(int i=0;i<nt;i++){
            int k = i;
            slegnpu_(rthk,rvs,rrho,rmax,t + i,cp+i,cout+k,uz,tz,
                    t1 + i,c1 + i, t2 + i,c2 + i,dcdbr,dcdhr,dcdrr,
                    dudbr,dudhr,dudrr,iflsph);
            _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dudbr,kvs+k*nz);
            _ParamConvert(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dudrr,krho+k*nz);
            // zero out kvp
            for(int iz = 0; iz < nz; iz ++) kvp[k * nz + iz] = 0.;
        } 
    }

    // free space
    delete[] rrho;delete[] rvp;delete[] rvs;delete[] rthk;
    delete[] ur; delete [] uz; delete [] tr; delete[] tz;
    delete [] dudrr; delete[] dudar; delete[] dudbr; delete[] dudhr;
    delete [] dcdrr; delete[] dcdar; delete[] dcdbr; delete[] dcdhr;
}


void compute_kernel_fd(float *vs,float *vp,float *rho,float *dep,int nz,
                    const double *t,int nt,double *cout,double *kvs,double *kvp,
                    double *krho,std::string &wavetype)
{
    // convert grid-based model to layer-based model
    const int sublayer = 3;
    int rmax = nz + (nz-1) * sublayer;
    float *rvp = new float [rmax],*rvs = new float [rmax];
    float *rthk = new float [rmax],*rrho = new float [rmax];
    float dv = 0.01;
    double c1[nt],c2[nt];

    // vp 
    for(int iter = 0; iter < 3; iter ++) {
        float *v;
        double *ker;
        if(iter == 0) {
            v = vp;
            ker = kvp;
        }
        else if (iter == 1) {
            v = vs;
            ker = kvs;
        }
        else {
            v = rho;
            ker = krho;
        }
        for(int iz = 0; iz < nz; iz ++) {
            float v0 = v[iz];
            v[iz] = v0 * ( 1 - 0.5  * dv);
            grid2LayerModel(dep,vp,vs,rho,nz,sublayer,rthk,rvp,rvs,rrho);
            surfdisp(rthk,rvp,rvs,rrho,rmax,t,c1,nt,wavetype,0,true,false);
            v[iz] = v0 * ( 1 + 0.5  * dv);
            grid2LayerModel(dep,vp,vs,rho,nz,sublayer,rthk,rvp,rvs,rrho);
            surfdisp(rthk,rvp,rvs,rrho,rmax,t,c2,nt,wavetype,0,true,false);

            for(int it = 0; it < nt; it ++) {
                ker[it * nz + iz] = (c2[it] - c1[it]) / (v0 * dv);
            }
            v[iz] = v0;
        }
        v = nullptr;
        ker  = nullptr;
    }

    delete[] rrho;delete[] rvp;delete[] rvs;delete[] rthk;
}


int main () {
    fvec z,vp,vs,rho;

    // read model
    int nz;
    FILE *fp = fopen("model.txt","r");
    fscanf(fp,"%d",&nz);
    z.resize(nz); vp.resize(nz); vs.resize(nz); rho.resize(nz);
    for(int i = 0; i < nz; i ++) {
        fscanf(fp,"%f%f%f%f",&z[i],&vs[i],&vp[i],&rho[i]);
    }
    fclose(fp);

    // period 
    int nt = 19;
    dvec t(nt),cg(nt);
    for(int it = 0; it < nt; it ++) t[it] = 4.0 + 2.0 * it;
    dmat2 kvp(nz,nt), kvs(nz,nt), krho(nz,nt);
    dmat2 kvp_fd(nz,nt), kvs_fd(nz,nt), krho_fd(nz,nt);

    // compute kernel 
    std::vector<std::string> wavetypes{"Rc","Rg","Lc","Lg"};
    for(int i = 0; i < 4; i ++) {
        printf("compute wavetype = %s ...\n",wavetypes[i].c_str());
        compute_kernel(vs.data(),vp.data(),rho.data(),z.data(),nz,t.data(),
                        nt,cg.data(),kvs.data(),kvp.data(),krho.data(),wavetypes[i]);
        compute_kernel_fd(vs.data(),vp.data(),rho.data(),z.data(),nz,t.data(),
                        nt,cg.data(),kvs_fd.data(),kvp_fd.data(),krho_fd.data(),wavetypes[i]);
        
        // save file
        for(int it = 0; it < nt; it ++) {
            std::string outname = "output/disper_" + wavetypes[i] + "_" +  std::to_string(it)  + ".txt";
            FILE *fp = fopen(outname.c_str(),"w");
            for(int iz = 0; iz < nz; iz ++) {
                fprintf(fp,"%f %g %g %g %g %g %g\n",z[iz],kvs(iz,it),kvs_fd(iz,it),krho(iz,it),krho_fd(iz,it),
                                            kvp(iz,it),kvp_fd(iz,it));
            }
            fclose(fp);
        }
    }
    return 0;
}