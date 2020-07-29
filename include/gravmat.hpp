#pragma once
#include<iostream>
#include<string>
#include"coo_matrix.hpp"
//=== Random observation system, all variables are in Observation-Centred Coordinate System
class OBSSphGraRandom{
    public:
    int   np;     //=== number of observation grids in x, y, directions ===
	float  z0;	//=== observation elevation
    float  *lon, *lat; //=== observation lon/lat

	float *Gr; //=== Observered gravtiy in radial direction.

    int  israd; //=== Flag for lon/lat (ISRAD=0) or lonrad/colatrad/r (ISRAD=1)====

    ~OBSSphGraRandom(){
        delete[] Gr; delete[] lon; delete[] lat;
    }

    void chancoor(int flag);
    void read_obs_data(std::string filename);
};

/*	This model is used to compute gravity effects*/
class MOD3DSphGra{
    public:
    int     nx, ny, nz;            		//=== number of inversion grids in x, y, z, directions ===
    float   x0, y0, z0;
    float   *lon, *lat, *dep;
	float   *density0, *density;
    int synflag;

	int	israd;				//===Flag for lon/lat/dep (israd=0) or lonrad/colatrad/r (israd=1)====

    ~MOD3DSphGra(){
        delete[] lon; delete[] lat; delete[] dep;
        delete[] density0; 
        if(synflag) delete[] density;
    }

    void chancoor(int flag);
    void read_model(std::string paramfile,std::string modinfile,std::string modtrue="nothing");
};

void gravmat(MOD3DSphGra &mod3dsphgra,OBSSphGraRandom &ObsSphGra,coo_matrix<float> &smat);
void gravmat_parallel(MOD3DSphGra &mod3dsphgra,OBSSphGraRandom &ObsSphGra,
                    coo_matrix<float> &smat);