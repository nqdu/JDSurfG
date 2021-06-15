#include"IOFunction.hpp"
#include<math.h>
#include"utils.hpp"
#include"tomography.hpp"
const int MAX_LEN = 500;
const double DEG2RAD = M_PI / 180.;

int read_receiver(FILE *fp,char *line,std::vector<float> &rcx,
                std::vector<float> &rcz,std::vector<float> &v)
{ 
    while(!feof(fp)){
        if(fgets(line,MAX_LEN*sizeof(char),fp)==NULL)
            break;
        if(line[0] == '#') break;
        float stalat,stalon,velvalue;
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

void __read_traveltime(std::string &surfdata,SurfTime &surf)
{
    char line[MAX_LEN];
    char dummy;
    float sta1_lat,sta1_lon;
    int wavetp,veltp,period;
    int dall=0;

    FILE *fp;
    if((fp=fopen(surfdata.c_str(),"r"))==NULL){
        std::cout <<"cannot open file" << std::endl;
        exit(0);
    }
    if(fgets(line,MAX_LEN*sizeof(char),fp)==NULL){
        std::cout <<"cannot read file" << std::endl;
        exit(0);
    }
    while(!feof(fp)){
        // extract source station information
        sscanf(line,"%c%f%f%d%d%d",&dummy,&sta1_lat,&sta1_lon,&period,&wavetp,&veltp);
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
        int nr = read_receiver(fp,line,rcx,rcz,v);

        // init station pair
        StationPair pair(wtp,dall,period-1,nr,sta1_lat,sta1_lon);
        for(int i=0;i<nr;i++){
            float dist;
            pair.rcx[i] = rcx[i];
            pair.rcz[i] = rcz[i];
            delsph(sta1_lat,sta1_lon,rcx[i],rcz[i],dist);
            pair.dist[i] = dist;
            pair.obstime[i] = dist / v[i];
            dall ++ ;
            //std::cout << pair.obstime[i] << std::endl;
        }
        surf.Pairs.push_back(pair);
    }

    // print data information on screen
    fclose(fp);
    std::cout<<"The number of traveltime measurements = "<<dall<<std::endl;
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
}

int __read_parfile(std::ifstream &infile,std::string &paramfile,
                    MOD3d &mod,SurfTime &surf,InverseParamsBase &param)
{
   // step1 : read paramfile
    std::string line; 
    std ::istringstream info;    

    // read model description
    skipread(infile,line,"%d%d%d",&mod.nx,&mod.ny,&mod.nz);
    skipread(infile,line,"%f%f",&mod.goxd,&mod.gozd);
    skipread(infile,line,"%f%f",&mod.dvxd,&mod.dvzd);

    // output some infomation to screen
    printf("Model Description:\n");
    printf("===================================\n");
    printf("model origin: latitude,longitude\n");
    printf("   %g   %g\n",mod.goxd,mod.gozd);
    printf("model grid spacing: dlat,dlon\n");
    printf("   %g   %g\n",mod.dvxd,mod.dvzd);
    printf("model dimension: nlat,nlon,nz\n");
    printf("%5d %5d %5d\n",mod.nx,mod.ny,mod.nz);

    // read forward parameters
    skipread(infile,line,"%d",&surf.sublayer);
    skipread(infile,line,"%d",&surf.num_threads);
    printf("\nForward/Adjoint Computation Parameters:\n");
    printf("===================================\n");
    printf("sublayer = %d\n",surf.sublayer);
    printf("Number of Threads Used = %d\n",surf.num_threads);

    // read inversion parameters
    int nthreads;
    skipread(infile,line,"%f%f",&param.smooth,&param.damp); // smooth and damp
    skipread(infile,line,"%f%f",&param.minvel,&param.maxvel); // min and max vel
    skipread(infile,line,"%d",&param.maxiter); // max iteration
    skipread(infile,line,"%d",&nthreads);
    printf("\nInversion Parameters:\n");
    printf("===================================\n");
    printf("smooth = %f, damp = %f\n",param.smooth,param.damp);
    printf("Min amd max velocity(km/s) = %f, %f\n",param.minvel,param.maxvel);
    printf("Max iterations = %d\n",param.maxiter);
    printf("Number of Threads Used in Solving Linear System = %d\n",nthreads);

    printf("\nDispersion Data:\n");
    printf("===================================\n");
    // read Rayleigh phase
    skipread(infile,line,"%d",&surf.kmaxRc);
    if(surf.kmaxRc > 0){
        skipread(infile,line);
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
    skipread(infile,line,"%d",&surf.kmaxRg);
    if(surf.kmaxRg > 0){
        skipread(infile,line);
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
    skipread(infile,line,"%d",&surf.kmaxLc);
    if(surf.kmaxLc > 0){
        skipread(infile,line);
        info.str(line);
        surf.tLc.resize(surf.kmaxLc);
        std::cout << "Love wave phase velocity used,periods:(s)" << std::endl;
        for(int i=0;i<surf.kmaxLc;i++){
            info >> surf.tLc(i);
            printf("%7.1f ",surf.tLc(i));
        }
        printf("\n");
    }

    // read love group
    skipread(infile,line,"%d",&surf.kmaxLg);
    if(surf.kmaxLg > 0){
        skipread(infile,line);
        info.str(line);
        surf.tLg.resize(surf.kmaxLg);
        std::cout << "Love wave group velocity used,periods:(s)" << std::endl;
        for(int i=0;i<surf.kmaxLg;i++){
            info >> surf.tLg(i);
            printf("%7.1f ",surf.tLg(i));
        }
        printf("\n");
    }

    // cache some temporary variables
    int kmaxRc,kmaxRg,kmaxLc,kmaxLg;
    kmaxRc = surf.kmaxRc;
    kmaxRg = surf.kmaxRg;
    kmaxLc = surf.kmaxLc;
    kmaxLg = surf.kmaxLg;
    
    // read synthetic test params
    skipread(infile,line,"%d",&param.ifsyn);

    // allocate data parameters
    int nx = mod.nx,ny = mod.ny, nz= mod.nz;
    surf.kmax = kmaxRc + kmaxRg + kmaxLc + kmaxLg;
    mod.lon.resize(ny); mod.lat.resize(nx);
    mod.dep.resize(nz);
    mod.vs.resize(nx,ny,nz);

    // renew lon and lat
    for(int i=0;i<ny;i++) mod.lon(i) = mod.gozd + i * mod.dvzd;
    for(int i=0;i<nx;i++) mod.lat(i) = mod.goxd - i * mod.dvxd;

    return nthreads;
}

void __read_InputModel(std::string &modfile,Eigen::Tensor<float,3> &vs,
                        float *dep,bool print_depth)
{
    // read open model file
    std::ifstream infile;
    infile.open(modfile);

    // temporary variables
    int nz = vs.dimension(2),nx = vs.dimension(0);
    int ny = vs.dimension(1);

    // read depth grids first
    for(int i=0;i<nz;i++){
        infile >> dep[i];
    }
    if(print_depth){
        printf("\nGrid points in depth direction:(km):\n");
        printf("===================================\n");
        //std::cout<<"grid points in depth direction:(km)" << std::endl;
        for(int i=0;i<nz;i++) printf("%7.2f ",dep[i]);
        std::cout<< std::endl;
    }

    // read initial S-wave model
    for(int k=0;k<nz;k++){
    for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++){
        infile >> vs(i,j,k);
    }}}
    infile.close();
}