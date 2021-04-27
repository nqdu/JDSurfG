#define EIGEN_DONT_PARALLELIZE
#include"tomography.hpp"
#include<fstream>
#include"utils.hpp"
#include"openmp.hpp"
using Eigen::Tensor;
using Eigen::VectorXf;
using Eigen::VectorXi;
#define DEG2RAD M_PI/180.

int read_receiver(FILE *fp,char *line,std::vector<float> &rcx,
                std::vector<float> &rcz,std::vector<float> &v)
{ 
    while(!feof(fp)){
        if(fgets(line,300*sizeof(char),fp)==NULL)
            break;
        if(line[0] == '#') break;
        float stalat,stalon,velvalue,dist1;
        sscanf(line,"%f%f%f",&stalat,&stalon,&velvalue);
        stalat=(90.0-stalat)* DEG2RAD;
        stalon=stalon* DEG2RAD;
        rcx.push_back(stalat);
        rcz.push_back(stalon);
        v.push_back(velvalue);
    } 
    int nr = rcx.size();

    return nr;   
}

// compute std of a vector
float dnrm2(float *a,int n)
{
    float mean = 0.0,s=0.0;
    for(int i=0;i<n;i++){
        mean += a[i];
        s += a[i] * a[i];
    }
    mean /= n;
    s /= n;

    float std = std - mean * mean;
    std = sqrt(std);

    return std;
}

void JointTomo::readdata(std::string paramfile,std::string modfile,std::string surfdata,
                 std::string gravdata,std::string gravmat,std::string refmod,
                 std::string modtrue)
{
   // step1 : read paramfile
    std::ifstream infile;
    std::string line;
    std ::istringstream info;    
    int nsrc;
    infile.open(paramfile);

    getline(infile,line);
    sscanf(line.c_str(),"%d%d%d",&mod.nx,&mod.ny,&mod.nz);

    getline(infile,line);
    sscanf(line.c_str(),"%f%f",&mod.goxd,&mod.gozd);

    getline(infile,line);
    sscanf(line.c_str(),"%f%f",&mod.dvxd,&mod.dvzd);

    getline(infile,line);
    sscanf(line.c_str(),"%f%f",&param.smooth,&param.damp);
 
    getline(infile,line);
    sscanf(line.c_str(),"%d",&surf.sublayer);

    getline(infile,line);
    sscanf(line.c_str(),"%f%f",&param.minvel,&param.maxvel);

    getline(infile,line);
    sscanf(line.c_str(),"%d",&param.maxiter);

    //getline(infile,line);
    //sscanf(line.c_str(),"%f",&param.spra);

    // output some infomation to screen
    printf("model origin: latitude,longitude\n");
    printf("%7.1f %7.1f\n",mod.goxd,mod.gozd);
    printf("model grid spacing: dlat,dlon\n");
    printf("   %g   %g\n",mod.dvxd,mod.dvzd);
    printf("model dimension: nlat,nlon,nz\n");
    printf("%5d %5d %5d\n",mod.nx,mod.ny,mod.nz);

    // read Rayleigh phase
    getline(infile,line);
    sscanf(line.c_str(),"%d",&surf.kmaxRc);
    if(surf.kmaxRc > 0){
        getline(infile,line);
        info.str(line);
        surf.tRc.resize(surf.kmaxRc);
        std::cout << "Rayleigh wave phase velocity used,periods:(s)" << std::endl;
        for(int i=0;i<surf.kmaxRc;i++){
            info >> surf.tRc(i);
            printf("%7.1f ",surf.tRc(i));
        }
        std::cout<<std::endl;
    }
    
    // read Rayleigh group
    getline(infile,line);
    sscanf(line.c_str(),"%d",&surf.kmaxRg);
    if(surf.kmaxRg > 0){
        getline(infile,line);
        info.str(line);
        surf.tRg.resize(surf.kmaxRg);
        std::cout << "Rayleigh wave group velocity used,periods:(s)" << std::endl;
        for(int i=0;i<surf.kmaxRg;i++){
            info >> surf.tRg(i);
            printf("%7.1f ",surf.tRg(i));
        }
        printf("\n");
    }

    // read love phase
    getline(infile,line);
    sscanf(line.c_str(),"%d",&surf.kmaxLc);
    if(surf.kmaxLc > 0){
        getline(infile,line);
        info.str(line);
        surf.tLc.resize(surf.kmaxLc);
        std::cout << "Love wave phase velocity used,periods:(s)" << std::endl;
        for(int i=0;i<surf.kmaxLc;i++){
            info >> surf.tLc(i);
            printf("%7.1f",surf.tLc(i));
        }
        printf("\n");
    }

    // read love group
    getline(infile,line);
    sscanf(line.c_str(),"%d",&surf.kmaxLg);
    if(surf.kmaxLg > 0){
        getline(infile,line);
        info.str(line);
        surf.tLg.resize(surf.kmaxLg);
        std::cout << "Love wave group velocity used,periods:(s)" << std::endl;
        for(int i=0;i<surf.kmaxLg;i++){
            info >> surf.tLg(i);
            printf("%7.1f",surf.tLg(i));
        }
        printf("\n");
    }

    int kmaxRc,kmaxRg,kmaxLc,kmaxLg;
    kmaxRc = surf.kmaxRc;
    kmaxRg = surf.kmaxRg;
    kmaxLc = surf.kmaxLc;
    kmaxLg = surf.kmaxLg;

    getline(infile,line);
    sscanf(line.c_str(),"%d",&param.ifsyn);

    getline(infile,line);
    sscanf(line.c_str(),"%f%f",&param.noiselevel,&param.noiselevel1);

    getline(infile,line);
    sscanf(line.c_str(),"%f",&param.p);

    getline(infile,line);
    sscanf(line.c_str(),"%f%f",&param.weight1,&param.weight2);

    // allocate data parameters
    int nx = mod.nx,ny = mod.ny, nz= mod.nz;
    surf.kmax = kmaxRc + kmaxRg + kmaxLc + kmaxLg;
    mod.dep.resize(nz);
    mod.lon.resize(ny); mod.lat.resize(nx);
    mod.vs.resize(nx,ny,nz);
    n  = (nx-2) * (ny -2) * (nz -1 );
    if(param.ifsyn) vstrue.resize(nx,ny,nz);
    // renew lon and lat
    for(int i=0;i<ny;i++) mod.lon(i) = mod.gozd + i * mod.dvzd;
    for(int i=0;i<nx;i++) mod.lat(i) = mod.goxd - i * mod.dvxd;

    infile.close();

    // step2 :read surface wave traveltime data
     // read surface wave traveltime data
    char line1[300];

    char dummy;
    float sta1_lat,sta1_lon,sta2_lat,sta2_lon,velvalue,dist1;
    int wavetp,veltp,period,maxnar;
    int dall=0;

    FILE *fp;
    if((fp=fopen(surfdata.c_str(),"r"))==NULL){
        std::cout <<"cannot open file" << std::endl;
        exit(0);
    }
    if(fgets(line1,300*sizeof(char),fp)==NULL){
        std::cout <<"cannot read file" << std::endl;
        exit(0);
    }
    while(!feof(fp)){

        // extract source station information
        sscanf(line1,"%c%f%f%d%d%d",&dummy,&sta1_lat,&sta1_lon,&period,&wavetp,&veltp);
        sta1_lat= (90.0-sta1_lat)*DEG2RAD;
        sta1_lon *= DEG2RAD;
        std::string wtp;
        if ( wavetp==2 && veltp==0 ) 
            wtp="Rc";
        else if ( wavetp==2 && veltp==1 ) 
            wtp="Rg";
        else if ( wavetp==1 && veltp==0 ) 
            wtp="Lc";
        else
            wtp="Lg";
        std::vector<float> rcx,rcz,v;
        int nr = read_receiver(fp,line1,rcx,rcz,v);

        // init station pair
        StationPair pair(wtp,dall,period-1,nr,sta1_lat,sta1_lon);
        for(int i=0;i<nr;i++){
            float dist;
            pair.rcx[i] = rcx[i];
            pair.rcz[i] = rcz[i];
            delsph(sta1_lat,sta1_lon,rcx[i],rcz[i],&dist);
            pair.dist[i] = dist;
            pair.obstime[i] = dist / v[i];
            dall ++ ;
            //std::cout << pair.obstime[i] << std::endl;
        }
        surf.Pairs.push_back(pair);
    }

    // print data information on screen
    fclose(fp);
    std::cout<<"The number of traveltime measurements is "<<dall<<std::endl;
    surf.num_data = dall;
    surf.obst.resize(dall);
    surf.sta_dist.resize(dall);
    int count = 0;
    for(int i=0;i<surf.Pairs.size();i++){
        for(int j=0;j<surf.Pairs[i].nr;j++){
            surf.obst(count) = surf.Pairs[i].obstime[j];
            surf.sta_dist(count) = surf.Pairs[i].dist[j];
            count += 1;
        }
    }

    // step3 : read gravity data and remove average value
    obsg.read_obs_data(gravdata);
    std::cout <<"The number of gravity measurements " <<obsg.np << std::endl;
    num_data = surf.num_data + obsg.np;

    // step4: read initial model
    infile.open(modfile);
    std::cout<<"grid points in depth direction:(km)" << std::endl;
    for(int i=0;i<mod.nz;i++){
        infile >> mod.dep(i);
        printf("%7.2f",mod.dep(i));
    }
    std::cout<< std::endl;

    // read initial S-wave model
    for(int k=0;k<mod.nz;k++){
    for(int j=0;j<mod.ny;j++){
    for(int i=0;i<mod.nx;i++){
        infile >> mod.vs(i,j,k);
    }}}
    infile.close();

    // step5 : read ture model if required
    if(param.ifsyn == 1){
        infile.open(modtrue);

        for(int k=0;k<mod.nz;k++){
        for(int j=0;j<mod.ny;j++){
        for(int i=0;i<mod.nx;i++){
            infile >> vstrue(i,j,k);
        }}}
        infile.close();
    }

    // step6: read refmodel
    modref = mod;
   if(refmod!="None"){
        std::cout<<"the gravity reference model is " + refmod <<std::endl;
        infile.open(refmod);
        for(int i=0;i<nz;i++) infile >>modref.vs(0,0,i);
        // read initial S-wave model
        for(int k=0;k<mod.nz;k++){
        for(int j=0;j<mod.ny;j++){
        for(int i=0;i<mod.nx;i++){
            infile >> modref.vs(i,j,k);
        }}}
        infile.close();
    }
    else{
        std::cout<<"no reference model is given, so the average of initial model is used" << std::endl;
        for(int k=0;k<mod.nz;k++){
            float mean = 0.0;
            for(int j=0;j<mod.ny-2;j++){
                for(int i=0;i<mod.nx-2;i++){
                    mean += mod.vs(i+1,j+1,k);
                }
            }
            mean /= (mod.nx-2) * (mod.ny - 2);
            for(int j=0;j<mod.ny;j++){
                for(int i=0;i<mod.nx;i++){
                    modref.vs(i,j,k) = mean;
                }
            }
        }
    }
    
    //step7 : read gravity matrix
    std::cout << "reading gravity matrix" << std::endl;
    int nar = 0;
    FILE *pin;
    if((pin=popen(("grep -v '#' " + gravmat + "| wc -l").c_str(), "r"))==NULL){
        std::cout << "cannot open file "<< gravmat << std::endl;
        exit(0);
    }
    int flag = fscanf(pin,"%d",&nar);
    pclose(pin);
    gmat.initialize(obsg.np,n,nar);
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
 * Assemble gravity and surfdisp matrix to a global one, and add weights
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
assemble(std::string basedir,Tensor<float,3> &vsf,csr_matrix<float> &smat,
         float weight1,float weight2)
{
    // get model and matrix dimensions
    int nx = mod.nx, ny = mod.ny, nz = mod.nz;
    int n = smat.cols();
    int ngrav = gmat.rows();
    int nsurf = surf.num_data;

    // assemble
    for(int r = 0;r<ngrav;r++){ // loop around all rows
        int start = gmat.indptr[r];
        int end = gmat.indptr[r + 1];
        int rwc = nsurf + r;
        smat.indptr[rwc + 1] = smat.indptr[rwc] + end - start;
    }
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
            empirical_relation(&vs,&vp,&rho);
            float tmp1,tmp2;
            empirical_deriv(vp,vs,&tmp1,&tmp2);

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
    int m = num_data;
    VectorXf res(m + n);

    // compute frechet kernel, and gravity anomaly
    std::string basedir = "kernelJ";
    VectorXi nonzeros = surf.FrechetKernel(mod,vsf,dsyn,basedir);
    modref.gravity(gmat,vsf,dg);
    float mean = dg.sum() / dg.size();
    dg.array() -= mean;

    // renew residues vector
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
    assemble(basedir,vsf,smat,sigma1,sigma2);

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
    LSMRDict<float> dict(itnlim,10,param.damp,param.smooth);
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
        vsf(i+1,j+1,k) = temp;
    }}}
}
