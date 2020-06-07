#include<fstream>
#include<sstream>
#include"defmod.hpp"
#include"empirical.hpp"
#include"spherical.hpp"
#include"const.hpp"

void OBSSphGraRandom :: chancoor(int flag)
{
    if( flag==0 && israd==1 ){          //change coordinates to lon/lat/dep;
        for(int i=0; i<np; i++ ){
            lon[i] = lon[i]/degrad2;
            lat[i] =  geocen2geogralat((double)(hpi2 - lat[i]) )/degrad2;
        }
        z0   = rearth2 - z0;
        israd = 0;
    }
    else if(flag==1 && israd==0) {     //change coordinates to lonrad/colatrad/r;
        for(int i=0; i<np; i++ ){
            lon[i] = lon[i]*degrad2;
            //lat[i] = hpi2 - lat[i]*degrad2;
            lat[i] = hpi2 - geogra2geocenlat( (double)(lat[i]*degrad2));
        }
        z0 = rearth2 - z0;
        israd = 1;
    }
	else{
        std::cout <<"don't need to change coordinates" << std::endl;
    }
}

void OBSSphGraRandom :: read_obs_data(std::string filename)
{
    std::ifstream infile;
    std::string line;
    np = 0;

    // get no. of lines
    infile.open(filename);
    while(!infile.eof()){
        getline(infile,line);
            np ++ ;
    }
    if(line == "") np --;
    infile.close();

    infile.open(filename);
    lon = new float [np];
    lat = new float [np];
    Gr = new float [np];
    z0 = 0.0;
    israd = 0;

    for(int i=0;i<np;i++){
        infile >> lon[i] >> lat[i] >> Gr[i];
    }

    infile.close();
}

