#include"JointTomo.hpp"
#include"delsph.hpp"
#include"calsurf.hpp"
#include"gaussian.hpp"
#include"empirical.hpp"
#include"const.hpp"
#include<fstream>
using namespace Eigen;

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
    for(int i=0;i<4;i++){
        getline(infile,line);
    }

    getline(infile,line);
    sscanf(line.c_str(),"%d%d%d",&mod.nx,&mod.ny,&mod.nz);

    getline(infile,line);
    sscanf(line.c_str(),"%f%f",&mod.goxd,&mod.gozd);

    getline(infile,line);
    sscanf(line.c_str(),"%f%f",&mod.dvxd,&mod.dvzd);

    getline(infile,line);
    sscanf(line.c_str(),"%d",&nsrc);

    getline(infile,line);
    sscanf(line.c_str(),"%f%f",&param.smooth,&param.damp);

    getline(infile,line);
    sscanf(line.c_str(),"%d",&surf.sublayer);

    getline(infile,line);
    sscanf(line.c_str(),"%f%f",&param.minvel,&param.maxvel);

    getline(infile,line);
    sscanf(line.c_str(),"%d",&param.maxiter);

    getline(infile,line);
    sscanf(line.c_str(),"%f",&param.spra);

    // output some infomation to screen
    printf("model origin: latitude,longitude\n");
    printf("%7.1f %7.1f\n",mod.goxd,mod.gozd);
    printf("model grid spacing: dlat,dlon\n");
    printf("   %g %g\n",mod.dvxd,mod.dvzd);
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
    sscanf(line.c_str(),"%f",&param.noiselevel);

    getline(infile,line);
    sscanf(line.c_str(),"%f",&param.p);

    getline(infile,line);
    sscanf(line.c_str(),"%f",&param.weight1);

    getline(infile,line);
    sscanf(line.c_str(),"%f",&param.weight2);

    // allocate data parameters
    int nx = mod.nx,ny = mod.ny, nz= mod.nz;
    surf.n = (nx-2) * (ny-2) * (nz -1);
    n = surf.n;
    surf.kmax = surf.kmaxRc + surf.kmaxRg + surf.kmaxLc + surf.kmaxLg;
    int kmax = surf.kmax;
    modref.nx = mod.nx; modref.ny = mod.ny;
    modref.nz = mod.nz;
    mod.dep.resize(nz);modref.dep.resize(nz);
    mod.lon.resize(ny); mod.lat.resize(nx);
    modref.lon.resize(ny); modref.lat.resize(nx);
    mod.vs.resize(nx,ny,nz);
    modref.vs.resize(nx,ny,nz);

    if(param.ifsyn) vstrue.resize(nx,ny,nz);
    surf.scxf.resize(nsrc,kmax);
    surf.sczf.resize(nsrc,kmax);
    surf.rcxf.resize(nsrc,nsrc,kmax);
    surf.rczf.resize(nsrc,nsrc,kmax);
    surf.periods.resize(nsrc,kmax);
    surf.wavetype.resize(nsrc,kmax);
    surf.nrc1.resize(nsrc,kmax);
    surf.nsrc1.resize(kmax);
    surf.igrt.resize(nsrc,kmax);

    // renew lon and lat
    for(int i=0;i<ny;i++) mod.lon(i) = mod.gozd + i * mod.dvzd;
    for(int i=0;i<nx;i++) mod.lat(i) = mod.goxd - i * mod.dvxd;

    infile.close();

    // step2 :read surface wave traveltime data
    char line1[300];
    float *obs = new float [nsrc * nsrc *kmax];
    float *dis = new float [nsrc* nsrc * kmax];

    int steprc=0,steps=0,knum=0,knumo;
    char dummy;
    float sta1_lat,sta1_lon,sta2_lat,sta2_lon,velvalue,dist1;
    int wavetp,veltp,period,maxnar;
    int dall=0;
    knumo=99999;

    FILE *fp;
    if((fp=fopen(surfdata.c_str(),"r"))==NULL){
        std::cout <<"cannot open file" << std::endl;
        exit(0);
    }
    while(!feof(fp)){
        // read one line
        if(fgets(line1,300*sizeof(char),fp)==NULL)
            break;

        // if read the source term
        if(line1[0]=='#'){
            sscanf(line1,"%c%f%f%d%d%d",&dummy,&sta1_lat,&sta1_lon,&period,&wavetp,&veltp);

            // change period index according to wavetype and velotype
            if ( wavetp==2 && veltp==0 ) 
                knum = period-1;
            else if ( wavetp==2 && veltp==1 ) 
                knum = period -1 +kmaxRc;
            else if ( wavetp==1 && veltp==0 ) 
                knum = period - 1 + kmaxRc + kmaxRg;
            else
                knum = period -1 + kmaxRc + kmaxRg + kmaxLc;
            
            if (knum!=knumo) // if get a new period, set number of source=0
                steps=0;
            steprc=0; // init number of receivers
            sta1_lat= (90.0-sta1_lat)*pi/180.0;
            sta1_lon *= pi/180.0;
            surf.scxf(steps,knum) = sta1_lat;
            surf.sczf(steps,knum) = sta1_lon;
            surf.periods(steps,knum) = period;
            surf.wavetype(steps,knum) = wavetp;
            surf.igrt(steps,knum) = veltp;
            surf.nsrc1(knum)=steps+1;
            knumo=knum;
            steps+=1;
        }
        else{
            sscanf(line1,"%f%f%f",&sta2_lat,&sta2_lon,&velvalue);
            sta2_lat=(90.0-sta2_lat)*pi/180.0;
            sta2_lon=sta2_lon*pi/180.0;
            surf.rcxf(steprc,steps-1,knum) = sta2_lat;
            surf.rczf(steprc,steps-1,knum) = sta2_lon;
            delsph(sta1_lat,sta1_lon,sta2_lat,sta2_lon,dis+dall);
            obs[dall] = dis[dall] / velvalue;
            dall++;
            surf.nrc1(steps-1,knum) = steprc + 1;
            steprc+=1;
        }
    }

    // print data information on screen
    fclose(fp);
    std::cout<<"The number of traveltime measurements is "<<dall<<std::endl;
    surf.num_data = dall;

    // allocate space for dist and obst
    surf.obst.resize(dall); surf.dist.resize(dall);
    for(int i=0;i<dall;i++){
        surf.obst(i) = obs[i];
        surf.dist(i) = dis[i];
    }
    delete[] dis; delete[] obs;

    // step3 : read gravity data
    obsg.read_obs_data(gravdata);
    std::cout <<"no. of gravity measurements " <<obsg.np << std::endl;
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
            }
        }
    }
    infile.close();

    // step5 : read ture model if required
    if(param.ifsyn){
        infile.open(modtrue);

        for(int k=0;k<mod.nz;k++){
            for(int j=0;j<mod.ny;j++){
                for(int i=0;i<mod.nx;i++){
                    infile >> vstrue(i,j,k);
                }
            }
        }
        infile.close();
    }

    // step6: read refmodel
   if(refmod!="None"){
        std::cout<<"the gravity reference model is " + refmod <<std::endl;
        infile.open(refmod);
        for(int i=0;i<nz;i++) infile >>modref.vs(0,0,i);
        // read initial S-wave model
        for(int k=0;k<mod.nz;k++){
            for(int j=0;j<mod.ny;j++){
                for(int i=0;i<mod.nx;i++){
                    infile >> modref.vs(i,j,k);
                }
            }
        }
        infile.close();
    }
    else{
        std::cout<<"no reference model is given, so the average of \
                    initial model is used" << std::endl;
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
    // read gravity matrix
    fp=fopen(gravmat.c_str(),"r");
    int nar = 0;
    std::cout << "reading gravity forward matrix ..."<<std::endl;
    while(!feof(fp)){
        // read one line
        if(fgets(line1,300*sizeof(char),fp)==NULL)
            break;
        nar ++;
    }
    fclose(fp);
    gmat.initialize(obsg.np,n,nar);
    fp=fopen(gravmat.c_str(),"r");
    for(int i=0;i<nar;i++){
        int t = fscanf(fp,"%d%d%f",gmat.rw+i,gmat.col+i,gmat.val+i);
    }
    fclose(fp);
}

void JointTomo ::forward(Tensor<float,3> &vs,VectorXf &dsyn,VectorXf &dg)
{   
    surf.forward(mod,vs,dsyn);
    modref.gravity(gmat,vs,dg);
}

void JointTomo:: checkerboard()
{
    VectorXf dsyn(surf.num_data),dg(obsg.np);
    forward(vstrue,dsyn,dg);

    // add noise
    for(int i=0;i<surf.num_data;i++){
        surf.obst(i) = dsyn(i)  +  param.noiselevel * gaussian();
    }

    for(int i=0;i<obsg.np;i++){
        obsg.Gr[i] = dg(i)  +  param.noiselevel1 * gaussian();
    }
}

void JointTomo:: inversion(Tensor<float,3> &vsf,VectorXf &dv,
                            VectorXf &dsyn,VectorXf &dg)
{
    int nx = mod.nx, ny = mod.ny;
    int nz = mod.nz;
    int m = num_data;
    VectorXf res(m + n);

    // initialize global sparse matrix
    int nonzeros =(int)(surf.num_data * n * param.spra) + gmat.nonzeros;
    coo_matrix<float>smat(m + n,n,nonzeros);

    // compute FrechetKernel
    int nar = surf.FrechetKernel(mod,vsf,dsyn,smat);
    modref.gravity(gmat,vsf,dg);
    // nar is no. of zeros of current smat, without regularization

    // compute residuals
    res.setZero();
    res.segment(0,surf.num_data) = surf.obst - dsyn;
    Map<VectorXf> Gr(obsg.Gr,obsg.np);
    res.segment(surf.num_data,obsg.np) = Gr - dg;

    // compute weights according to Julia(2000)
    float sigma1,sigma2;
    if(param.ifsyn == 1){
        sigma1 = 1.0;
        sigma2= 1.0;
    }
    else if(param.ifsyn ==2){
        float sigma1 = dnrm2(res.data(),m);
        float sigma2 = dnrm2(res.data() + m,obsg.np);
        sigma1 /= sigma2;
        sigma2 = 1.0;
    }
    else{
        sigma1 = param.weight1;
        sigma2 = param.weight2;
    }
    sigma1 = sqrt(param.p /surf.num_data) / sigma1;
    sigma2 = sqrt((1-param.p)/obsg.np) / sigma2;

    // copy gravity matrix to joint matrix
    // empirical relation is used to convert dg/drho to dg/dvs
    // convert change of density to change of S wave velocity
    // some other empirical relation maybe better
    // here we use:
    //! rho = 1.6612*vp-0.4721*vp**2 + &
    //            0.0671*vp**3 - 0.0043*vp**4 + 0.000106*vp**5
    // vp = 0.9409 + 2.0947*vs - 0.8206*vsz**2+ 0.2683*vs**3 - 0.0251*vs**4
    // then G*drho = g0 could be changed to G'*dvs = g0 
    for(int p=0;p<gmat.nonzeros;p++){
        // compute index corresponding to current column
	    float vp,vs,rho;
        int col = gmat.col[p];
        int k = col /((nx-2) * (ny-2));
        int c = col %((nx-2) * (ny-2));
        int i = c%(nx-2);
        int j = c/(nx-2);

        // compute drho/dvs, store it in tmp1 * tmp2
        vs = vsf(i,j,k);
        empirical_relation(&vs,&vp,&rho);
        float tmp1,tmp2;
        empirical_deriv(vp,vs,&tmp1,&tmp2);

        // append gmat to smat
        smat.val[nar + p] = gmat.val[p] * tmp1 * tmp2;
        smat.col[nar + p] = col;
        smat.rw[nar + p] = surf.num_data + gmat.rw[p];
    }

    // add weights
    for(int p=0;p<nar;p++){
        smat.val[p] *= sigma1;
    }
    for(int p=0;p<gmat.nonzeros;p++){
        smat.val[p+nar] *= sigma2;
    } 
    for(int i=0;i<surf.num_data;i++){
        res(i) *= sigma1;
    }
    for(int i=0;i<obsg.np;i++){
        res(i+surf.num_data) *= sigma2;
    }
    nar += gmat.nonzeros;

    // add regularization terms
    int count = 0;
    float smooth = param.smooth;
    for(int k=0;k<nz-1;k++){
    for(int j=0;j<ny-2;j++){
    for(int i=0;i<nx-2;i++){
        if( i==0 || i==nx-3 || j==0 || j==ny-3 || k==0 || k==nz-2){
            // absorb large variations at boundary points
            if(nar + 1 > smat.nonzeros){
                std::cout << "please increase sparse ratio!" << std::endl;
                exit(0);
            }
            smat.col[nar] = k * (ny -2) * (nx -2) + j * (nx -2) + i;
            smat.val[nar] = 2.0 * smooth;
            smat.rw[nar] = count + m;
            nar ++;
            count ++ ; 
        }
        else{
            if(nar  + 7 > smat.nonzeros){
                std::cout << "please increase sparse ratio!" << std::endl;
                exit(0);
            }
            int rwc = count + m;  // current row
            int clc = k * (ny -2) * (nx -2) + j * (nx -2) + i;// current column
            smat.val[nar] = 6.0 * smooth;
            smat.col[nar] = clc;
            smat.rw[nar] = rwc;

            // x direction
            smat.val[nar + 1] = -smooth;
            smat.rw[nar + 1] = rwc;
            smat.col[nar + 1] = clc - 1;
            smat.val[nar + 2] = -smooth;
            smat.rw[nar + 2] = rwc;
            smat.col[nar + 2] = clc + 1;

            // y direction
            smat.val[nar + 3] = -smooth;
            smat.rw[nar + 3] = rwc;
            smat.col[nar + 3] = clc - (nx - 2);
            smat.val[nar + 4] = -smooth;
            smat.rw[nar + 4] = rwc;
            smat.col[nar + 4] = clc + (nx - 2);

            // z direction
            smat.val[nar + 5] = -smooth;
            smat.rw[nar + 5] = rwc;
            smat.col[nar + 5] = clc - (nx - 2) * (ny - 2);
            smat.val[nar + 6] = -smooth;
            smat.rw[nar + 6] = rwc;
            smat.col[nar + 6] = clc + (nx - 2) * (ny - 2);

            nar += 7;
            count++;
        }
    }}}

    // renew sparse matrix meta-informations
    //smat.m = m + count;
    smat.nonzeros = nar;
    smat.cpp2fortran();

   // solve equations by lsmr
    std::cout <<"solving linear systems by LSMR ..." << std::endl;
    int itnlim = n * 2;
    LSMRDict<float> dict(itnlim,10,param.damp,param.smooth);
    smat.LsmrSolver(res.data(),dv.data(),dict);
    std::cout << "max negative and positive perturbation: " \
                << dv.minCoeff() <<" " << dv.maxCoeff()\
                <<std::endl;

    // renew vsf and tackle large variations
    TensorMap<Tensor<float,3>>dx(dv.data(),nx-2,ny-2,nz-1);
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
        vsf(i+1,j+1,k) = temp;
    }}}
}