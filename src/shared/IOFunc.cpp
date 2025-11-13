#include <iostream>
#include <dirent.h>
#include <sys/stat.h>
#include <fstream>

#include "numerical.hpp"

/**
 * @brief Create a directory if it doesn not exist
 * 
 * @param dirname 
 */
void 
create_directory(const char *dirname)
{
    if(opendir(dirname) == NULL){
        mkdir(dirname,0777);
    } 
}

/**
 * @brief print progress bar to the screen
 * @param percentage percentage of current progress
*/
void print_progressbar(float percentage) {
    // Ensure percentage is within [0.0, 100.0]
    percentage = (percentage < 0.0) ? 0.0 : (percentage > 100.0) ? 100.0 : percentage;

    // Determine the number of characters to represent the progress bar
    int numChars = (int)(percentage / 2.0);

    // Print the progress bar
    printf("[");
    for (int i = 0; i < 50; ++i) {
        if (i < numChars) {
            printf("=");
        } else if (i == numChars) {
            printf(">");
        } else {
            printf(" ");
        }
    }
    printf("] %.1f%%\r", percentage);  // Use carriage return to overwrite the line

    // Flush the output to ensure immediate display
    fflush(stdout);
}

/**
 * @brief Read a velocity model from a file
 * 
 * @param modfile The path to the model file
 * @param vs Reference to a 3D float matrix to store the velocity model
 * @param depth Reference to a float vector to store depth values
 * @param goxd Reference to a float to store the origin latitude
 * @param gozd Reference to a float to store the origin longitude
 * @param dvxd Reference to a float to store the latitude grid spacing
 * @param dvzd Reference to a float to store the longitude grid spacing
 * @param print_info Flag to indicate whether to print information to the screen
 */
void read_velocity_model(const std::string &modfile,
                         fmat3 &vs,
                         fvec &depth,
                         fvec &lon,
                         fvec &lat,
                         bool print_info)
{
    int nx,ny,nz;
    float goxd,gozd,dvxd,dvzd;
    std::ifstream infile; infile.open(modfile);
    if(!infile.is_open()) {
        printf("cannot open %s\n",modfile.c_str());
        exit(1);
    }
    std::string line;
    getline(infile,line);
    sscanf(line.c_str(),"%d%d%d",&nx,&ny,&nz);
    getline(infile,line);
    sscanf(line.c_str(),"%f%f",&goxd,&gozd);
    getline(infile,line);
    sscanf(line.c_str(),"%f%f",&dvxd,&dvzd);

    // set lat and lon vectors
    lon.resize(ny); lat.resize(nx);
    for(int i = 0; i < nx; i ++){
        lat[i] = goxd - i * dvxd;
    }
    for(int i = 0; i < ny; i ++){
        lon[i] = gozd + i * dvzd;
    }

    // output some infomation to screen
    if(print_info) {
        printf("\nModel Description:\n");
        printf("===================================\n");
        printf("model origin: latitude,longitude\n");
        printf("   %g   %g\n",goxd,gozd);
        printf("model grid spacing: dlat,dlon\n");
        printf("   %g   %g\n",dvxd,dvzd);
        printf("model dimension: nlat,nlon,nz\n");
        printf("%5d %5d %5d\n",nx,ny,nz); 
    }

    // allocate space
    vs.resize(nx,ny,nz);
    depth.resize(nz);
    
    // read depth
    if(print_info) {
        printf("Grid points in depth direction:(km):\n");
    }
    getline(infile,line);
    size_t len = line.size();
    char tmp[len + 10];
    strcpy(tmp,line.c_str());
    char *starp = tmp,*endp = NULL;
    for(int i = 0; i < nz; i ++) {
        depth[i] = std::strtof(starp,&endp);
        starp = endp;
        if(print_info) {
            printf("%7.2f ",depth[i]);
        }
    }
    if(print_info) {
        printf("\n\n");
    }

    // read model
    for(int k=0;k<nz;k++){
    for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++){
        infile >> vs(i,j,k);
    }}}
    infile.close();
}

/**
 * @brief Interpolate topography data onto the model grid
 * 
 * @param topofile topofiles
 * @param lat latitude
 * @param lon longitude 
 * @param depth local 1-D depth model to be updated
 */
void interpolate_topo(const std::string &topofile,
                     const fvec &lat,
                     const fvec &lon,
                     fmat3 &depth)
{
    int nx = depth.dimension(0);
    int ny = depth.dimension(1);
    int nz = depth.dimension(2);

    // read topo
    FILE *fp = fopen(topofile.c_str(),"r");
    if(fp == NULL) {
        printf("cannot open %s\n",topofile.c_str());
        exit(1);
    
    }
    int nx_topo, ny_topo;
    float lon0,lat0, dlon, dlat;
    fscanf(fp,"%d%d",&nx_topo,&ny_topo);
    fscanf(fp,"%f%f",&lat0,&lon0);
    fscanf(fp,"%f%f",&dlat,&dlon);
    fmat2 topo(ny_topo,nx_topo);
    for(int j=0;j<ny_topo;j++){
    for(int i=0;i<nx_topo;i++){
        fscanf(fp,"%f",&topo(i,j));
    }}
    fclose(fp);

    // loop over all grid points to interpolate
    for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++){
        // bilinear interpolation
        float lati = lat[j];
        float loni = lon[i];

        // Handle edge cases: clamp indices to valid range
        int i0 = std::max(0, std::min((int)((lat0 - lati) / dlat), nx_topo - 2));
        int j0 = std::max(0, std::min((int)((loni - lon0) / dlon), ny_topo - 2));
        int i1 = i0 + 1;
        int j1 = j0 + 1;
        float q11 = topo(i0,j0);
        float q21 = topo(i1,j0);
        float q12 = topo(i0,j1);
        float q22 = topo(i1,j1);
        float x1 = lon0 + j0 * dlon;
        float x2 = lon0 + j1 * dlon;
        float y1 = lat0 - i0 * dlat;
        float y2 = lat0 - i1 * dlat;
        float topo_ij = (q11 * (x2 - loni) * (y2 - lati) +
                        q21 * (loni - x1) * (y2 - lati) +
                        q12 * (x2 - loni) * (lati - y1) +
                        q22 * (loni - x1) * (lati - y1)) / ((x2 - x1) * (y2 - y1));
        
        depth(i,j,0) = -topo_ij;

        // make sure the thickness is non-negative
        float thk = depth(i,j,1) - depth(i,j,0);
        if(thk < 0.0) {
            printf("Warning: negative thickness at lat=%g, lon=%g, set to %g km\n",
                   lati,loni,thk);
            thk = 0.01;
            printf("reset topography from %g to %g km\n",depth(i,j,0),depth(i,j,1)-thk);
        }
    }}
}