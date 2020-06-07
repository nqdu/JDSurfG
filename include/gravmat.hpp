#pragma once
#include<iostream>
#include<string>
#include"defmod.hpp"

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
                    coo_matrix<float> &smat,int nthreads = 4);