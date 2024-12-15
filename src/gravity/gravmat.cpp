/**
 * This script is the core module of adaptive Gauss-Legendre gravity modeling
 * Written by Zhiwei Li at Dec. 2011
 * Modified by Nanqiao Du by adding parallel module and coordinate format matrix
*/
#include "gravity/gravity_module.hpp"
#include <cmath>
#include "shared/parallel.hpp"
#include "gravity/transform.hpp"
#include "gravity/gauleg.hpp"
#include "shared/IOFunc.hpp"
#include "shared/csr_matrix.hpp"

const double  degrad2 = M_PI/180.0;
const double  hpi2 =  M_PI*0.5;
const float  rearth2=6371.0;
const float G = 6.667E-6;
/*
	Change NMAX and NSCALE bigger for more accurate results,
	but bigger NMAX and NSCALE means much more computation.
*/
//===minimum and maximum nodes for GLQ
//===the ratio of distance/GLQ grid spacing
//===the minmum scale for GLQ, in km
const int NMIN=5,NMAX=256;
const int NSCALE = 5;
const double MINDISKM = 0.5;
const float maxdis = 100.0;// maximum distance between 2 points 
const float ftol = 1.0e-6;

void OBSSphGraRandom :: chancoor(int flag)
{
    if( flag==0 && israd==1 ){          //change coordinates to lon/lat/dep;
        for(int i=0; i<np; i++ ){
            lon[i] = lon[i]/degrad2;
            lat[i] =  geocen2geogralat((double)(hpi2 - lat[i]) )/degrad2;
            z0[i] = rearth2 - z0[i];
        }
        //z0   = rearth2 - z0;
        israd = 0;
    }
    else if(flag==1 && israd==0) {     //change coordinates to lonrad/colatrad/r;
        for(int i=0; i<np; i++ ){
            lon[i] = lon[i]*degrad2;
            //lat[i] = hpi2 - lat[i]*degrad2;
            lat[i] = hpi2 - geogra2geocenlat( (double)(lat[i]*degrad2));
            z0[i] = rearth2 - z0[i];
        }
        //z0 = rearth2 - z0;
        israd = 1;
    }
	else{
        std::cout <<"don't need to change coordinates" << std::endl;
    }
}

void OBSSphGraRandom :: read_obs_data(const std::string &filename)
{
    // open file
    std::ifstream infile;
    infile.open(filename);
    if(!infile.is_open()) {
        printf("cannot open %s\n",filename.c_str());
        exit(1);
    }

    // get # of lines
    std::string line;
    np = 0;
    while(std::getline(infile,line)) {
        np += 1;
    }
    infile.close();

    // allocate space
    lon = new float [np];
    lat = new float [np];
    Gr = new float [np];
    z0 = new float[np]();
    israd = 0;

    // open again
    infile.open(filename);
    for(int i = 0; i < np; i ++) {
        infile >> lon[i] >> lat[i] >> Gr[i];
    }
    infile.close();

    printf("\nGravity Data: %d points\n",np);
}

void MOD3DSphGra :: chancoor(int flag)
{
	int	i;
	if( flag==0 && israd==1 ){//change coordinates to lon/lat/dep;
		for(i=0; i<nx; i++ )
            lon[i] = lon[i]/degrad2;

        for(i=0; i<ny; i++ )
			lat[i] = geocen2geogralat((double)(hpi2 - lat[i]) )/degrad2;

        for(i=0; i<nz; i++)
            dep[i] = rearth2 - dep[i];

		x0	 = x0/degrad2;
		y0   = geocen2geogralat((double)(hpi2 - y0) )/degrad2;

        z0   = rearth2 - z0;

        israd = 0;
	}
	else if( flag==1 && israd==0) {//change coordinates to lonrad/colatrad/r;
		for(i=0; i<nx; i++ )
            lon[i] = lon[i]*degrad2;

        for(i=0; i<ny; i++ )
			lat[i] =   hpi2 - geogra2geocenlat( (double)(lat[i]*degrad2) );


        for(i=0; i<nz; i++)
            dep[i] = rearth2 - dep[i];

		x0   = x0*degrad2;
		y0   = hpi2 - geogra2geocenlat( (double)(y0*degrad2) );
        z0   = rearth2 - z0;

		israd = 1;
	}
	else{
        std::cout <<"don't need to change coordinates" << std::endl;
	}
}

void MOD3DSphGra :: 
read_model(const std::string &modfile)
{
    // open file
    std::ifstream infile;
    infile.open(modfile);
    if(!infile.is_open()) {
        printf("cannot open %s\n",modfile.c_str());
        exit(1);
    }

    int i;
    float ulx,uly,dx,dy;
    std::string line;

    // read parameter file
    std::getline(infile,line);
    sscanf(line.c_str(),"%d%d%d",&ny,&nx,&nz);
    nx = nx-2;
    ny = ny-2;
    nz = nz-1;
    lon = new float[nx];
    lat = new float[ny];
    dep = new float[nz];
    std::getline(infile,line);
    sscanf(line.c_str(),"%f%f",&uly,&ulx);
    std::getline(infile,line);
    sscanf(line.c_str(),"%f%f",&dy,&dx);

    // remove boundaries
    ulx = ulx+dx;
    uly = uly-dy;
    x0 = ulx;
    y0 = uly;
    for(i = 0; i < nx; i ++) lon[i] = ulx + i * dx;
    for(i = 0; i < ny; i ++) lat[i] = uly - i * dy;

    israd = 0;

    // print model information
    printf("Model Description:\n");
    printf("===================================\n");
    printf("model origin: latitude,longitude\n");
    printf("   %g   %g\n",uly,ulx);
    printf("model grid spacing: dlat,dlon\n");
    printf("   %g   %g\n",dy,dx);
    printf("model dimension: nlat,nlon,nz\n");
    printf("%5d %5d %5d\n",ny,nx,nz);

    // read depth of velomodel
    std::getline(infile,line);
    char *tmp = new char[line.size() + 1];
    strcpy(tmp,line.c_str());
    char *starp = tmp, *endp = NULL;
    for(int i = 0; i < nz; i ++) {
        dep[i] = std::strtof(starp,&endp);
        starp = endp;
    }
    delete[] tmp;

    // close file
    infile.close();
}


/*
    Given position in lonrad/colatrad/r, get the density or the anomally density
*/
void getdensityBlock(MOD3DSphGra &mod3dsphgra, int px, int py, int pz, double *density )
{
    int     nx, ny, nz;
    long    N, nxy;
                                                                                                                                                                                
    nx = mod3dsphgra.nx;
    ny = mod3dsphgra.ny;
    nz = mod3dsphgra.nz;
    nxy = nx*ny;
                                                                                                                                                                                   
    if( px<0 || px>=nx || py<0 || py>=ny || pz<0 || pz>=nz ){
            fprintf(stderr,"point out of model. px, py, pz=%d %d %d\n", px, py, pz );
            fprintf(stderr,"model:nx, ny, nz = %d %d %d\n", nx, ny, nz );
            exit(0);
    }
                                                                                                                                                                                
    N = pz * nxy + px * ny + py;
    *density = (double)(mod3dsphgra.density[N]-mod3dsphgra.density0[N]);
}

/*
	for Radial Gravity, from GravityEffects
*/
void GravityR(double ro,double lonrado,double colatrado,
            double rs,double lonrads,double colatrads, 
            double density,double *Gr)
{       
        double  K, costmp, Ros, Ros3;
        
        //====change to International Unit kg/m3
        density = density*1000.0;
        ro      = ro*1000.0;
        rs      = rs*1000.0;
        
        K = G*density*rs*rs*sin(colatrads);
        costmp = cos(colatrado)*cos(colatrads) + sin(colatrado)*sin(colatrads)*cos(lonrado-lonrads );
        Ros = sqrt(ro*ro + rs*rs - 2.0*ro*rs*costmp );
        
        Ros3 = Ros*Ros*Ros;
        
         //NOTE: According the Observation coordinate, the Downward/Northward/Eastward of the observation station
        // on or above the Earth's Surface,
        // Chnage the Polarity for: Gr, Gcolat, Grlon, Gcolatlon, Glonr, Gloncolat

        (*Gr) = (-1.0) *K*( -ro + rs*costmp )/Ros3;
}

/*
    Given one observation point outside the model, compute the Gravity matrix from all the cells of model.
    Using the adaptive scheme according the distance from the source to the observation point.
*/
void gravmat_one_row(MOD3DSphGra &mod3dsphgra, double ro, double lonrado,
                     double colatrado,int *col,float *value,int &nar)
{
    double  gr;
    double  rs1, lonrads1, colatrads1, rs2, lonrads2, colatrads2;
    double  rs, lonrads, colatrads;
    double  dlonrad, dcolatrad, dr;
    double  alonrad, acolatrad, ar;
    double   weight, wgl,Gr;

    int     nx, ny, nz;
    int     i, j, k;
    int     ii, jj, kk;
    int n;

    /*
            parameters for Gauss-Legendre nodes number
    */
    double   mindis;
    float   dx, dy ,dz;
    int     xn, yn, zn;
    double  x1 = -1.0, x2 = 1.0;
    double  xw[NMAX+1], yw[NMAX+1], zw[NMAX+1];   //NOTICE:       xw, yw, zw, xp, yp, zp begin from 1 to xn, yn, zn
    double  xp[NMAX+1], yp[NMAX+1], zp[NMAX+1];   //NOTICE:       NOT from 0.
    k = NMAX+1 ;
    for(i=0; i<k; i++){
        xw[i] = 0.0; yw[i] = 0.0; zw[i] = 0.0;
        xp[i] = 0.0; yp[i] = 0.0; zp[i] = 0.0;
    }

    nx = mod3dsphgra.nx;
    ny = mod3dsphgra.ny;
    nz = mod3dsphgra.nz;

    nar = 0;
    for( k=(nz-1); k>=0; k-- ){//===from bottom to top==
    for( i=0; i<nx; i++ ){
    for( j=0; j<ny; j++ ){
        Gr = 0.0;
        n = k*nx*ny+i*ny+j;
        //get the beginning and ending points for each direction.
        if( k == (nz-1) ){
            rs1 = mod3dsphgra.dep[k];
            rs2 = 0.5*(mod3dsphgra.dep[k-1]+mod3dsphgra.dep[k]);
        }
        else if ( k == 0 ){
            rs1 = 0.5*(mod3dsphgra.dep[k]+mod3dsphgra.dep[k+1]);
            rs2 = mod3dsphgra.dep[k];
        }
        else{
            rs1 = 0.5*(mod3dsphgra.dep[k]+mod3dsphgra.dep[k+1]);
            rs2 = 0.5*(mod3dsphgra.dep[k-1]+mod3dsphgra.dep[k]);
        }
                                                                                                                                                                                                
        if( i == 0 ){
            lonrads1 = mod3dsphgra.lon[i];
            lonrads2 = 0.5*(mod3dsphgra.lon[i+1]+mod3dsphgra.lon[i]);
        }
        else if( i == (nx-1) ){
            lonrads1 = 0.5* (mod3dsphgra.lon[i]+mod3dsphgra.lon[i-1]);
            lonrads2 = mod3dsphgra.lon[i];
        }
        else{
            lonrads1 = 0.5*(mod3dsphgra.lon[i]+mod3dsphgra.lon[i-1]);
            lonrads2 = 0.5*(mod3dsphgra.lon[i+1]+mod3dsphgra.lon[i]);
        }
                                                                                                                                                                                                
        if( j == 0 ){
            colatrads1 = mod3dsphgra.lat[j];
            colatrads2 = 0.5*(mod3dsphgra.lat[j+1]+mod3dsphgra.lat[j]);
        }
        else if ( j == (ny-1) ){
            colatrads1 = 0.5*(mod3dsphgra.lat[j]+mod3dsphgra.lat[j-1]);
            colatrads2 = mod3dsphgra.lat[j];
        }
        else{
            colatrads1 = 0.5*(mod3dsphgra.lat[j]+mod3dsphgra.lat[j-1]);
            colatrads2 = 0.5*(mod3dsphgra.lat[j+1]+mod3dsphgra.lat[j]);
        }
                                                                                                                                                                                                
        //getdensityBlock( mod3dsphgra, i, j, k, &density );      //===get the density at the grid i/j/k block

        //====change to International Unit kg/m3
        //weight        = 1000.0*((rs2-rs1)*(lonrads2-lonrads1)*(colatrads2-colatrads1))/8.0;
        weight  = 125.0*((rs2-rs1)*(lonrads2-lonrads1)*(colatrads2-colatrads1));

        dlonrad = lonrads2 - lonrads1;
        dcolatrad = colatrads2 - colatrads1;
        dr = rs2 - rs1;
        alonrad = lonrads2 + lonrads1;
        acolatrad = colatrads2 + colatrads1;
        ar = rs2 + rs1;

        //====calculate the distance from source to observation point,
        //====determine the parameters for Gauss-Legendre quadrature integration
        dx = fabs( dlonrad*6371.0*sin(colatrads1) );
        dy = fabs( dcolatrad *6371.0 );
        dz = fabs( dr );

        cal_minDis( mindis, colatrado, lonrado, ro, colatrads1, lonrads1, colatrads2, lonrads2, rs2 );
        if( mindis < MINDISKM )  mindis = MINDISKM;
        if(mindis<maxdis){

            xn = (int) ( NSCALE/(mindis/dx) );//====set the scale of each small block is 1/10 of the distance to the observation point
            yn = (int) ( NSCALE/(mindis/dy) );
            zn = (int) ( NSCALE/(mindis/dz) );
            if( xn < NMIN ) xn = NMIN;
            if( yn < NMIN ) yn = NMIN;
            if( zn < NMIN ) zn = NMIN;
            if( xn > NMAX ) xn = NMAX;
            if( yn > NMAX ) yn = NMAX;
            if( zn > NMAX ) zn = NMAX;

            gauleg( x1, x2, xp+1, xw+1,  xn);
            gauleg( x1, x2, yp+1, yw+1,  yn);
            gauleg( x1, x2, zp+1, zw+1,  zn);

            //===do Gauss-Legendre quadrature integration===
            for( kk=1; kk<=zn; kk++ ){
            for( jj=1; jj<=yn; jj++ ){
            for( ii=1; ii<=xn; ii++ ){
                lonrads         = (xp[ii]*dlonrad + alonrad)/2.0;
                colatrads       = (yp[jj]*dcolatrad + acolatrad)/2.0;
                rs              = (zp[kk]*dr + ar )/2.0;
                wgl             = weight*xw[ii]*yw[jj]*zw[kk];

                //getdensity( mod3dsphgra, lonrads, colatrads, rs, &density );
                GravityR( ro, lonrado, colatrado, rs, lonrads, colatrads, 1.0, &gr );
                Gr += wgl*gr;
            }}}
            if(Gr > ftol){
                value[nar] = Gr;
                col[nar] = n;
                nar += 1;
            }
        }
    }}}

}

/*
    Given Random Observation System. Only Compute matrix of gravity in Radial direction.
    openmp is used.
*/
void generate_gravmat(MOD3DSphGra &mod3dsphgra,OBSSphGraRandom &ObsSphGra,
                      const std::string &outfile)
{
    int np,m;
    np  = ObsSphGra.np;
    m = mod3dsphgra.nx * mod3dsphgra.ny * mod3dsphgra.nz;

    // get num threads
    int nprocs{};
    #pragma omp parallel
    {
        nprocs = omp_get_num_threads();
    }

    #pragma omp parallel for shared(ObsSphGra,mod3dsphgra)
    for(int myrank = 0; myrank < nprocs; myrank ++) {
        // open tempory file to save gravmat 
        std::string filename = outfile + "." + std::to_string(myrank);
        std::ofstream fp(filename,std::ios::binary);

        // write rows and cols of this matrix
        fp.write((char*)&np,sizeof(int));
        fp.write((char*)&m,sizeof(int));

        int starid,endid;
        int col[m];
        float data[m];
        allocate_tasks(np,nprocs,myrank,starid,endid);
        for(int i = starid; i <= endid; i ++) {
            double lonrado   = ObsSphGra.lon[i];
            double colatrado = ObsSphGra.lat[i];
            double ro = ObsSphGra.z0[i];
            int nar;
            gravmat_one_row(mod3dsphgra,ro,lonrado,colatrado,col,data,nar);
            
            if(myrank == 0) {
                int size = endid - starid + 1;
                float per = (i-starid+1.) / size * 100;
                print_progressbar(per);
                // printf("computing %d of %d point: lon=%f lat=%f\n",
                //         i-starid+1,size,lonrado/degrad2,90.-colatrado/degrad2);
            }

            // write matrix to temporary file
            fp.write((char*)&i,sizeof(int));
            fp.write((char*)&nar,sizeof(int));
            fp.write((char*)col,sizeof(int) * nar);
            fp.write((char*)data,sizeof(float) * nar);
        }

        // close temporary file
        fp.close();
    }

    // merge files to one big file
    merge_csr_files(nprocs,outfile);
}

