#include"tomography.hpp"
#include"utils.hpp"
#include"coo_matrix.hpp"
#include<fstream>
#include"openmp.hpp"
using Eigen::Tensor;
using Eigen::VectorXf;

int read_receiver(FILE *fp,char *line,std::vector<float> &rcx,
                std::vector<float> &rcz,std::vector<float> &v)
{
    while(!feof(fp)){
        if(fgets(line,300*sizeof(char),fp)==NULL)
            break;
        if(line[0] == '#') break;
        float stalat,stalon,velvalue,dist1;
        sscanf(line,"%f%f%f",&stalat,&stalon,&velvalue);
        stalat=(90.0-stalat)*pi/180.0;
        stalon=stalon*pi/180.0;
        rcx.push_back(stalat);
        rcz.push_back(stalon);
        v.push_back(velvalue);
    } 
    int nr = rcx.size();

    return nr;   
}

int SurfTomo:: readdata(std::string paramfile,std::string datafile,
                        std::string modfile,std::string modtrue )
{
 // read paramfile
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
    sscanf(line.c_str(),"%f",&param.noiselevel);

    // allocate data parameters
    int nx = mod.nx,ny = mod.ny, nz= mod.nz;
    surf.kmax = kmaxRc + kmaxRg + kmaxLc + kmaxLg;
    mod.dep.resize(nz);
    mod.lon.resize(ny); mod.lat.resize(nx);
    mod.vs.resize(nx,ny,nz);
    unknowns = (nx-2) * (ny -2) * (nz -1 );
    if(param.ifsyn) vstrue.resize(nx,ny,nz);

    // renew lon and lat
    for(int i=0;i<ny;i++) mod.lon(i) = mod.gozd + i * mod.dvzd;
    for(int i=0;i<nx;i++) mod.lat(i) = mod.goxd - i * mod.dvxd;

    infile.close();

    // read surface wave traveltime data
    char line1[300];

    char dummy;
    float sta1_lat,sta1_lon,sta2_lat,sta2_lon,velvalue,dist1;
    int wavetp,veltp,period,maxnar;
    int dall=0;

    FILE *fp;
    if((fp=fopen(datafile.c_str(),"r"))==NULL){
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
        sta1_lat= (90.0-sta1_lat)*pi/180.0;
        sta1_lon *= pi/180.0;
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
    std::cout<<"The number of measurements is "<<dall<<std::endl;
    surf.num_data = dall;
    num_data = dall;
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

    // initialize sparse matrix
    //int nonzeros = (int)(dall * n * param.spra + n);
    //smat.initialize(surf.num_data+n,n,nonzeros);

    // read initial model
    infile.open(modfile);

    // read depth grids first
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

    // read ture model if required
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

    return 1;
}

void SurfTomo:: forward(Eigen::Tensor<float,3> &vs,Eigen::VectorXf &data)
{
    surf.TravelTime(mod,vs,data);
}

int SurfTomo::FrechetKernel(Tensor<float,3> &vs,VectorXf &data,std::string save_dir)
{
    int nar = surf.FrechetKernel(mod,vs,data,save_dir);

    return nar;
}

void SurfTomo :: checkerboard()
{
    VectorXf dsyn(surf.num_data);
    forward(vstrue,dsyn);
    // add noise
    for(int i=0;i<surf.num_data;i++){
        dsyn(i) *= (1.0 + param.noiselevel * gaussian());
    }
    surf.obst = dsyn;
}

void SurfTomo::inversion(Tensor<float,3> &vsf,VectorXf &dsyn)
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
    int nar = FrechetKernel(vsf,dsyn,basedir);
    //std::cout << nar << std::endl;
    res.segment(0,m) = surf.obst - dsyn;

    // initialize matrix
    coo_matrix<float> smat(m+n,n,nar + n*7);
    smat.setZeros();
    int count = 0;
    std::cout << "Assembling derivative Matrix ..." << std::endl;
    for(int i=0;i<nthread;i++){
        std::string filename = basedir + "/" +  std::to_string(i) + ".txt";
        count = smat.read(filename,count);
    }
    
    // add regularization terms
    count = 0;
    float weight = param.smooth;
    for(int k=0;k<nz-1;k++){
    for(int j=0;j<ny-2;j++){
    for(int i=0;i<nx-2;i++){
        if( i==0 || i==nx-3 || j==0 || j==ny-3 || k==0 || k==nz-2){
            // and more restrictions to boundary points
            if(nar + 1 > smat.nonzeros){
                std::cout << "please increase sparse ratio!" << std::endl;
                exit(0);
            }
            smat.col[nar] = k * (ny -2) * (nx -2) + j * (nx -2) + i;
            smat.val[nar] = 2.0 * weight;
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
            smat.val[nar] = 6.0 * weight;
            smat.col[nar] = clc;
            smat.rw[nar] = rwc;

            // x direction
            smat.val[nar + 1] = -weight;
            smat.rw[nar + 1] = rwc;
            smat.col[nar + 1] = clc - 1;
            smat.val[nar + 2] = -weight;
            smat.rw[nar + 2] = rwc;
            smat.col[nar + 2] = clc + 1;

            // y direction
            smat.val[nar + 3] = -weight;
            smat.rw[nar + 3] = rwc;
            smat.col[nar + 3] = clc - (nx - 2);
            smat.val[nar + 4] = -weight;
            smat.rw[nar + 4] = rwc;
            smat.col[nar + 4] = clc + (nx - 2);

            // z direction
            smat.val[nar + 5] = -weight;
            smat.rw[nar + 5] = rwc;
            smat.col[nar + 5] = clc - (nx - 2) * (ny - 2);
            smat.val[nar + 6] = -weight;
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
    LSMRDict<float> dict(itnlim,10,param.damp,weight);
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
        if(temp > maxvel && maxvel > 0.0) temp = maxvel;
        if(temp < minvel && minvel > 0.0) temp = minvel;
        vsf(i+1,j+1,k) = temp;
    }}}
}