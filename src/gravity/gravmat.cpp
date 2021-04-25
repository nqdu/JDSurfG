/**
 * This script is the core module of adaptive Gauss-Legendre gravity modeling
 * Written by Zhiwei Li at Dec. 2011
 * Modified by Nanqiao Du by adding parallel module and coordinate format matrix
*/
#include"utils.hpp"
#include"openmp.hpp"
#include"gravmat.hpp"
#include<fstream>
#include<stdlib.h>
#include<math.h>
#define  degrad2 (M_PI/180.0)
#define  hpi2    (M_PI*0.5)
#define  rearth2 6371.0
#define G 6.667E-6
/*
	Change NMAX and NSCALE bigger for more accurate results,
	but bigger NMAX and NSCALE means much more computation.
*/
//===minimum and maximum nodes for GLQ
//===the ratio of distance/GLQ grid spacing
//===the minmum scale for GLQ, in km
#define	NMIN      5
#define NMAX      256
#define NSCALE    5
#define MINDISKM  0.5
#define maxdis 100.0 // maximum distance between 2 points 
#define ftol 1.0e-6

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

    // get lines of this file
    FILE *pin;
    if((pin=popen(("wc -l " + filename).c_str(), "r"))==NULL){
        std::cout << "cannot open file "<< filename << std::endl;
        exit(0);
    }
    int flag = fscanf(pin,"%d",&np);
    pclose(pin);

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

void MOD3DSphGra :: read_model(std::string paramfile,std::string modfile)
{
    std::ifstream fp;
    int i,j,k;
    int n,kmax;
    float ulx,uly,dx,dy;
    float v,vtrue;
    std::string line;

    // read parameter file
    fp.open(paramfile);
    getline(fp,line);
    sscanf(line.c_str(),"%d%d%d",&ny,&nx,&nz);
    nx=nx-2;
    ny=ny-2;
    nz=nz-1;
    n = nx * ny * nz;
    lon = new float[nx];
    lat = new float[ny];
    dep = new float[nz];
    getline(fp,line);
    sscanf(line.c_str(),"%f%f",&uly,&ulx);
    getline(fp,line);
    sscanf(line.c_str(),"%f%f",&dy,&dx);
    ulx=ulx+dx;
    uly=uly-dy;
    x0=ulx;
    y0=uly;

    for(i=0;i<nx;i++)
        lon[i]=ulx+i*dx;
    for(i=0;i<ny;i++)
        lat[i]=uly-i*dy;

    //read synflag
    for(i=0;i<5;i++)
        getline(fp,line);
    sscanf(line.c_str(),"%d",&kmax); //kmaxRc
    if(kmax>0){
        getline(fp,line);
        getline(fp,line);
    }
    else{
        getline(fp,line);
    }

    sscanf(line.c_str(),"%d",&kmax); //kmaxRg
     if(kmax>0){
        getline(fp,line);
        getline(fp,line);
    }
    else{
        getline(fp,line);
    } 

    sscanf(line.c_str(),"%d",&kmax); //kmaxLc
     if(kmax>0){
        getline(fp,line);
        getline(fp,line);
    }
    else{
        getline(fp,line);
    }     

    sscanf(line.c_str(),"%d",&kmax); //kmaxLg
     if(kmax>0){
        getline(fp,line);
        getline(fp,line);
    }
    else{
        getline(fp,line);
    }   
    fp.close();
    israd = 0;

    // read depth of velomodel
    fp.open(modfile);
    for(i=0;i<nz;i++){
        fp >> dep[i];
    }
    /*
    fp >> v;

    // allocate space for density
    density0 = new float [n];
    
    // velocity model has shape (dep,lon,lat)
    float vp;
    for(k = 0;k < nz+1;k++){
    for(i = 0;i < nx+2;i++){
    for(j=0;j<ny+2;j++){
        n = k*nx * ny + (i-1) * ny + j -1;
        fp >> v;
        if(i == 0||j == 0||k == nz||i == nx+1||j == ny+1)
            continue;
        empirical_relation(&v,&vp,density0+n);
    }}}
    israd = 0;*/
    fp.close();
}


/*
    Given position in lonrad/colatrad/r, get the density or the anomally density
*/
void getdensityBlock(MOD3DSphGra &mod3dsphgra, int px, int py, int pz, double *density )
{
    int     nx, ny, nz;
    int     i, j, k;
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
                     double colatrado,int *col,float *value,int *nar)
{
    double  gr;
    double  rs1, lonrads1, colatrads1, rs2, lonrads2, colatrads2;
    double  rs, lonrads, colatrads;
    double  dlonrad, dcolatrad, dr;
    double  alonrad, acolatrad, ar;
    double   density;
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
                value[*nar] = Gr;
                col[*nar] = n;
                (*nar) += 1;
            }
        }
    }}}

}

/*
    Given Random Observation System. Only Compute matrix of gravity in Radial direction.
*/
void gravmat(MOD3DSphGra &mod3dsphgra,OBSSphGraRandom &ObsSphGra,csr_matrix<float> &smat)
{
    double   ro, lonrado, colatrado;
    int     i, j;
    int     np,nar,nar1;
    int     *col;
    int nx ,ny,nz;

    // extract parameters
    np  = ObsSphGra.np;
    ro = ObsSphGra.z0;
    nar = 0;
    nar1 = 0;
    nx = mod3dsphgra.nx;
    ny = mod3dsphgra.ny;
    nz = mod3dsphgra.nz;

    // arrays to store all the corresponding columns
    col = new int [(int)(0.1 * np * nx * ny * nz)]();

    for( i = 0; i < np; i++ ){
        lonrado   = ObsSphGra.lon[i];
        colatrado = ObsSphGra.lat[i];
        gravmat_one_row(mod3dsphgra,ro,lonrado,colatrado,col,smat.data,&nar1);
        smat.indptr[i+1] = smat.indptr[i] + nar1 - nar;
        printf("computing %d-th point: lon=%f lat=%f\n",
              i+1,lonrado/degrad2,90-colatrado/degrad2);
        nar = nar1;
    }

    for(i = 0;i < nar;i++){
        smat.indices[i] = col[i];
    }

    delete[] col;
    smat.nonzeros = nar;
}

/*
    Given Random Observation System. Only Compute matrix of gravity in Radial direction.
    openmpi is used.
*/
void gravmat_parallel(MOD3DSphGra &mod3dsphgra,OBSSphGraRandom &ObsSphGra,
                    csr_matrix<float> &smat)
{
    double   ro;
    int     i, j;
    int     np,nar,m;
    int nx,ny,nz;

    np  = ObsSphGra.np;
    ro = ObsSphGra.z0;
    m = mod3dsphgra.nx * mod3dsphgra.ny * mod3dsphgra.nz;

    // allocate sparse matrix for every data
    int *col,*nars;
    float *data;
    int ncol = (int)(m * 0.2);

    col = new int [np * ncol]();
    data = new float[np * ncol]();
    nars = new int [np]();

    omp_set_num_threads(nthreads);
    #pragma omp parallel for shared(ObsSphGra,mod3dsphgra)
    for(i=0;i<np;i++){
        double lonrado   = ObsSphGra.lon[i];
        double colatrado = ObsSphGra.lat[i];
        gravmat_one_row(mod3dsphgra,ro,lonrado,colatrado,
                        col+i*ncol,data+i*ncol,nars+i);
        printf("computing %d-th point: lon=%f lat=%f\n",
              i+1,lonrado/degrad2,90-colatrado/degrad2);
    }
    nar = 0; 

    for(i=0;i<np;i++){
        smat.indptr[i+1] = smat.indptr[i] + nars[i];
        for(j=0;j<nars[i];j++){
            smat.indices[nar] = col[i*ncol+j];
            smat.data[nar] = data[i*ncol+j];
            nar += 1;
        }
    }
    smat.nonzeros = nar; 

    delete[] col;
    delete[] nars;
    delete[] data;

}

#undef	NMIN
#undef NMAX
#undef NSCALE
#undef MINDISKM
#undef maxdis
#undef ftol
#undef  degrad2
#undef  hpi2
#undef  rearth2
#undef  G

