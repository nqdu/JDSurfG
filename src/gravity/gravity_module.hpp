#ifndef JDSURFG_GRAVITY_GRAVITY_MODULE_H_
#define JDSURFG_GRAVITY_GRAVITY_MODULE_H_

#include <iostream>
#include <string>

//=== Random observation system, all variables are in Observation-Centred Coordinate System
class OBSSphGraRandom{
    public:
    int   np;     //=== number of observation grids in x, y, directions ===
	//float  z0;	//=== observation elevation
    float  *lon, *lat,*z0; //=== observation lon/lat/elevation

	float *Gr; //=== Observered gravtiy in radial direction.

    int  israd; //=== Flag for lon/lat (ISRAD=0) or lonrad/colatrad/r (ISRAD=1)====

    ~OBSSphGraRandom(){
        delete[] Gr; delete[] lon; 
        delete[] lat; delete[] z0;
    }

    void chancoor(int flag);
    void read_obs_data(const std::string &filename);
};

/*	This model is used to compute gravity effects*/
class MOD3DSphGra{
    public:
    int     nx, ny, nz;            		//=== number of inversion grids in x, y, z, directions ===
    float   x0, y0, z0;
    float   *lon, *lat, *dep;
	float   *density0, *density;
    //int synflag;

	int	israd;				//===Flag for lon/lat/dep (israd=0) or lonrad/colatrad/r (israd=1)====

    ~MOD3DSphGra(){
        delete[] lon; delete[] lat; delete[] dep;
        //delete[] density0; 
        //if(synflag) delete[] density;
    }

    void chancoor(int flag);
    void read_model(const std::string &modinfile);
};

void generate_gravmat(MOD3DSphGra &mod3dsphgra,OBSSphGraRandom &ObsSphGra,
                      const std::string &outfile);

#endif // end JDSURFG_GRAVITY_GRAVITY_MODULE_H_
