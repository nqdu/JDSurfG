#include"surftomo.hpp"
#include<fstream>
#include"delsph.hpp"
#include"const.hpp"
using namespace Eigen;

void SurfTomo ::forward(Tensor<float,3> &vs,VectorXf &dsyn)
{
    surf.forward(mod,vs,dsyn);
}

void SurfTomo::checkerboard()
{
    VectorXf dsyn(surf.num_data);
    surf.checkerboard(mod,vstrue,dsyn,param.noiselevel);
    surf.obst = dsyn;
}

void SurfTomo :: inversion(Tensor<float,3> &vs,VectorXf &dv,VectorXf &dsyn)
{
    int nar = (int)(num_data * n * param.spra);
    coo_matrix<float> smat(num_data+ n,n,nar);
    surf.inversion(mod,vs,smat,dv,dsyn,param.damp,
                param.smooth,param.minvel,param.maxvel);
}

void SurfTomo ::readdata(std::string paramfile,std::string datafile,
                std::string modfile,std::string modtrue )
{
   // read paramfile
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
    surf.n = (nx-2) * (ny-2) * (nz -1);
    n = surf.n;
    surf.kmax = surf.kmaxRc + surf.kmaxRg + surf.kmaxLc + surf.kmaxLg;
    int kmax = surf.kmax;
    mod.dep.resize(nz);
    mod.lon.resize(ny); mod.lat.resize(nx);
    mod.vs.resize(nx,ny,nz);
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

    // read surface wave traveltime data
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
    if((fp=fopen(datafile.c_str(),"r"))==NULL){
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
    std::cout<<"The number of measurements is "<<dall<<std::endl;
    surf.num_data = dall;
    num_data = dall;

    // allocate space for dist and obst
    surf.obst.resize(dall); surf.dist.resize(dall);
    for(int i=0;i<dall;i++){
        surf.obst(i) = obs[i];
        surf.dist(i) = dis[i];
    }
    delete[] dis; delete[] obs;

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
}