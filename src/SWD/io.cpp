#include "SWD/SurfaceWave.hpp"
#include "numerical.hpp"
#include "shared/spherical.hpp"
#include <fstream>
const int MAX_LEN = 8192;

static int 
read_receiver(FILE *fp,char *line,std::vector<float> &rcx,
              std::vector<float> &rcz,std::vector<float> &v)
{ 
    const float DEG2RAD = M_PI / 180.;
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

void SurfTime::  
read_swd_data(const std::string &datafile)
{
    const float DEG2RAD = M_PI / 180.;
    char line[MAX_LEN];
    char dummy;
    float sta1_lat,sta1_lon;
    int wavetp,veltp,period,mode;
    int dall=0;

    FILE *fp;
    if((fp=fopen(datafile.c_str(),"r"))==NULL){
        printf("cannot open file %s\n",datafile.c_str());
        exit(1);
    }
    printf("Data Description:\n");
    printf("===================================\n");

    // read period vectors
    this -> kmax = 0;

    // Rc
    assert(fgets(line,MAX_LEN*sizeof(char),fp)!=NULL);
    char *starp = line,*endp = NULL;
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
        assert(fgets(line,MAX_LEN*sizeof(char),fp)!=NULL);
        starp = line; endp = NULL;
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
    assert(fgets(line,MAX_LEN*sizeof(char),fp)!=NULL);
    starp = line; endp = NULL;
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
        assert(fgets(line,MAX_LEN*sizeof(char),fp)!=NULL);
        starp = line; endp = NULL;
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
    assert(fgets(line,MAX_LEN*sizeof(char),fp)!=NULL);
    starp = line; endp = NULL;
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
        assert(fgets(line,MAX_LEN*sizeof(char),fp)!=NULL);
        starp = line; endp = NULL;
        for(int it = 0; it < nt; it ++){
            this ->tLc(it,im) = std::strtof(starp,&endp);
            starp = endp;
        }
    }

    // Lg
    assert(fgets(line,MAX_LEN*sizeof(char),fp)!=NULL);
    starp = line; endp = NULL;
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
        assert(fgets(line,MAX_LEN*sizeof(char),fp)!=NULL);
        starp = line; endp = NULL;
        for(int it = 0; it < nt; it ++){
            this ->tLg(it,im) = std::strtof(starp,&endp);
            starp = endp;
        }
    }

    if(fgets(line,MAX_LEN*sizeof(char),fp)==NULL){
        printf("cannot read file %s\n",datafile.c_str());
        exit(1);
    }
    while(!feof(fp)){
        // extract source station information
        sscanf(line,"%c%f%f%d%d%d%d",&dummy,&sta1_lat,&sta1_lon,&mode,&period,&wavetp,&veltp);
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
    fclose(fp);
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