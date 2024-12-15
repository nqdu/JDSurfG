#include "SWD/SurfaceWave.hpp"
#include "numerical.hpp"
#include "shared/spherical.hpp"
#include <fstream>

static int 
read_receiver(std::ifstream &infile,std::string &line,std::vector<float> &rcx,
              std::vector<float> &rcz,std::vector<float> &v)
{ 
    const float DEG2RAD = M_PI / 180.;
    while(!infile.eof()) {
        if(!std::getline(infile,line)) {
            break;
        }
        if(line[0] == '#') break;
        float stalat,stalon,velvalue;
        sscanf(line.c_str(),"%f%f%f",&stalat,&stalon,&velvalue);
        stalat=(90.0-stalat)* DEG2RAD;
        stalon=stalon* DEG2RAD;
        rcx.push_back(stalat);
        rcz.push_back(stalon);
        v.push_back(velvalue);
    } 
    int nr = rcx.size();

    return nr;   
}

void SurfTime::  
read_swd_data(const std::string &datafile)
{
    const float DEG2RAD = M_PI / 180.;
    std::string line;
    char dummy;
    float sta1_lat,sta1_lon;
    int wavetp,veltp,period,mode;
    int dall = 0;

    // tmp char
    std::vector<char> tmp;
    tmp.reserve(8192);

    // open file
    std::ifstream infile;
    infile.open(datafile);
    if(!infile.is_open()){
        printf("cannot open file %s\n",datafile.c_str());
        exit(1);
    }
    printf("Data Description:\n");
    printf("===================================\n");

    // read period vectors
    this -> kmax = 0;

    // Rc
    std::getline(infile,line);
    tmp.resize(line.size() + 1); memcpy(tmp.data(),line.c_str(),line.size());
    char *starp = tmp.data(), *endp = NULL;
    int nmode = std::strtof(starp,&endp); starp = endp;
    this ->kmaxRc.resize(nmode);
    for(int i = 0; i < nmode; i ++){
        this ->kmaxRc[i] = std::strtof(starp,&endp); starp = endp;
        this ->kmax += this ->kmaxRc[i];
    }
    int nt = 0;
    if(nmode > 0) nt = this ->kmaxRc.maxCoeff();
    this ->tRc.resize(nt,nmode); this ->tRc.setConstant(-1);
    for(int im = 0; im < nmode; im ++) {
        std::getline(infile,line);
        tmp.resize(line.size() + 1); memcpy(tmp.data(),line.c_str(),line.size());
        starp = tmp.data(); endp = NULL;
        for(int it = 0; it < nt; it ++){
            this ->tRc(it,im) = std::strtof(starp,&endp);
            starp = endp;
        }
    }
    if(nmode > 0) {
        printf("Rayleigh wave phase velocity used,periods:(s)\n");
        for(int im = 0; im < nmode; im ++) {
            printf("mode %d: ",im);
            for(int it = 0; it < nt; it ++){
                printf("%7.1f ",this ->tRc(it,im));
            }
            printf("\n");
        }
        printf("\n");
    }

    // Rg
    std::getline(infile,line);
    tmp.resize(line.size() + 1); memcpy(tmp.data(),line.c_str(),line.size());
    starp = tmp.data(); endp = NULL;
    nmode = std::strtof(starp,&endp); starp = endp;
    this ->kmaxRg.resize(nmode);
    for(int i = 0; i < nmode; i ++){
        this ->kmaxRg[i] = std::strtof(starp,&endp); starp = endp;
        this ->kmax += this ->kmaxRg[i];
    }
    nt = 0;
    if(nmode > 0) nt = this ->kmaxRg.maxCoeff();
    this ->tRg.resize(nt,nmode); this ->tRg.setConstant(-1);
    for(int im = 0; im < nmode; im ++) {
        std::getline(infile,line);
        tmp.resize(line.size() + 1); memcpy(tmp.data(),line.c_str(),line.size());
        starp = tmp.data(); endp = NULL;
        for(int it = 0; it < nt; it ++){
            this ->tRg(it,im) = std::strtof(starp,&endp);
            starp = endp;
        }
    }
    if(nmode > 0) {
        printf("Rayleigh wave group velocity used,periods:(s)\n");
        for(int im = 0; im < nmode; im ++) {
            printf("mode %d: ",im);
            for(int it = 0; it < nt; it ++){
                printf("%7.1f ",this ->tRg(it,im));
            }
            printf("\n");
        }
        printf("\n");
    }

    // Lc
    std::getline(infile,line);
    tmp.resize(line.size() + 1); memcpy(tmp.data(),line.c_str(),line.size());
    starp = tmp.data(); endp = NULL;
    nmode = std::strtof(starp,&endp); starp = endp;
    this ->kmaxLc.resize(nmode);
    for(int i = 0; i < nmode; i ++){
        this ->kmaxLc[i] = std::strtof(starp,&endp); starp = endp;
        this ->kmax += this ->kmaxLc[i];
    }
    nt = 0;
    if(nmode > 0) nt = this ->kmaxLc.maxCoeff();
    this ->tLc.resize(nt,nmode); this ->tLc.setConstant(-1);
    for(int im = 0; im < nmode; im ++) {
        std::getline(infile,line);
        tmp.resize(line.size() + 1); memcpy(tmp.data(),line.c_str(),line.size());
        starp = tmp.data(); endp = NULL;
        for(int it = 0; it < nt; it ++){
            this ->tLc(it,im) = std::strtof(starp,&endp);
            starp = endp;
        }
    }
    if(nmode > 0) {
        printf("Love wave phase velocity used,periods:(s)\n");
        for(int im = 0; im < nmode; im ++) {
            printf("mode %d: ",im);
            for(int it = 0; it < nt; it ++){
                printf("%7.1f ",this ->tLc(it,im));
            }
            printf("\n");
        }
        printf("\n");
    }

    // Lg
    std::getline(infile,line);
    tmp.resize(line.size() + 1); memcpy(tmp.data(),line.c_str(),line.size());
    starp = tmp.data(); endp = NULL;
    nmode = std::strtof(starp,&endp); starp = endp;
    this ->kmaxLg.resize(nmode);
    for(int i = 0; i < nmode; i ++){
        this ->kmaxLg[i] = std::strtof(starp,&endp); starp = endp;
        this ->kmax += this ->kmaxLg[i];
    }
    nt = 0;
    if(nmode > 0) nt = this ->kmaxLg.maxCoeff();
    this ->tLg.resize(nt,nmode); this ->tLg.setConstant(-1.);
    for(int im = 0; im < nmode; im ++) {
        std::getline(infile,line);
        tmp.resize(line.size() + 1); memcpy(tmp.data(),line.c_str(),line.size());
        starp = tmp.data(); endp = NULL;
        for(int it = 0; it < nt; it ++){
            this ->tLg(it,im) = std::strtof(starp,&endp);
            starp = endp;
        }
    }
    if(nmode > 0) {
        printf("Love wave group velocity used,periods:(s)\n");
        for(int im = 0; im < nmode; im ++) {
            printf("mode %d: ",im);
            for(int it = 0; it < nt; it ++){
                printf("%7.1f ",this ->tLg(it,im));
            }
            printf("\n");
        }
        printf("\n");
    }

    // read swd data 
    std::getline(infile,line);
    while(!infile.eof()) {
        // extract source station information
        sscanf(line.c_str(),"%c%f%f%d%d%d%d",&dummy,&sta1_lat,&sta1_lon,&mode,&period,&wavetp,&veltp);
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
        int nr = read_receiver(infile,line,rcx,rcz,v);

        // init station pair
        StationPair pair(wtp,dall,mode,period-1,nr,sta1_lat,sta1_lon);
        for(int i = 0;i < nr;i ++){
            float dist;
            pair.rcx[i] = rcx[i];
            pair.rcz[i] = rcz[i];
            dist = delsph(sta1_lat,sta1_lon,rcx[i],rcz[i]);
            pair.dist[i] = dist;
            pair.obstime[i] = dist / v[i];
            dall ++ ;
            //std::cout << pair.obstime[i] << std::endl;
        }
        this ->Pairs.push_back(pair);
    }

    // print data information on screen
    infile.close();
    printf("The number of traveltime measurements = %d\n",dall);
    this -> obst.resize(dall);
    this -> sta_dist.resize(dall);
    int count = 0;
    for(int i=0;i<this ->Pairs.size();i++){
        for(int j=0;j<this ->Pairs[i].nr;j++){
            this ->obst(count) = this ->Pairs[i].obstime[j];
            this ->sta_dist(count) = this ->Pairs[i].dist[j];
            count += 1;
        }
    }
}