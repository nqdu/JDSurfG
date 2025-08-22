#include "shared/csr_matrix.hpp"
#include "numerical.hpp"
#include "SWD/empirical.hpp"

#include <iostream>
#include <fstream>
#include <string>

/**
 * @brief compute gravity from a given model
 * 
 * @param vs vs model, shape(nx,ny,nz)
 * @param vsref reference vs model, shape(nx,ny,nz)
 * @param gmat gravity matrix, smat type
 * @param dsyn synthetic gravity data
 */
static void 
compute_gravity(const fmat3 &vs,const fmat3 &vsref,const csr_matrix &gmat,
                fvec &dsyn)
{
    int nx = vs.dimension(0), ny = vs.dimension(1), nz = vs.dimension(2);
    int n = (nx - 2) * (ny - 2) * (nz - 1);
    fvec drho(n); drho.setZero();
    dsyn.resize(gmat.rows()); dsyn.setConstant(0);

    // loop to compute drho
    int ic = 0;
    for(int k=0;k<nz-1;k++){
    for(int j=0;j<ny-2;j++){
    for(int i=0;i<nx-2;i++){
        float a,rho;
        empirical_relation(vs(i+1,j+1,k),a,rho);
        drho[ic] = rho;
        empirical_relation(vsref(i+1,j+1,k),a,rho);
        drho[ic] -= rho;
        ic += 1;
    }}}

    // compute gravity
    gmat.aprod(1,drho.data(),dsyn.data());
}

int main(int argc, char **argv){
    if(argc != 7) {
        printf("please run syngrav (reference model) (target model) obsgrav " \
                "gravmat remove_mean(0 or 1) outfile\n");
        printf("example: ./syngrav MOD.ref MOD.target obsgrav.dat gravmat.bin 1 out.txt\n");
        exit(1);
    }

    // get input args
    std::string refmod = std::string(argv[1]);
    std::string truemod = std::string(argv[2]);
    std::string gravfile = std::string(argv[3]);
    std::string gravmat = std::string(argv[4]);
    int remove_avg = std::stoi(argv[5]);
    std::string outfile = std::string(argv[6]);

    // variables needed
    fmat3 vsref,vs;
    int nx,ny,nz;

    // read reference model
    printf("reading ref/input model ...\n");
    std::ifstream infile; infile.open(refmod);
    std::string line;
    if(infile.fail()) {
        printf("cannot open %s\n",refmod.data());
        exit(1);
    }
    getline(infile,line);
    sscanf(line.data(),"%d%d%d",&nx,&ny,&nz);
    vs.resize(nx,ny,nz); vsref.resize(nx,ny,nz);
    for(int i = 0; i < 3; i ++)  getline(infile,line);

    for(int k=0;k<nz;k++){
    for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++){
        infile >> vsref(i,j,k);
    }}}
    infile.close();

    // read true model
    infile.open(truemod);
    if(infile.fail()) {
        printf("cannot open %s\n",truemod.data());
        exit(1);
    }
    for(int i = 0; i < 4; i ++) getline(infile,line);
    for(int k=0;k<nz;k++){
    for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++){
        infile >> vs(i,j,k);
    }}}
    infile.close();

    // read gravity matrix
    printf("reading gravity matrix ...\n");
    csr_matrix gmat;
    gmat.read_binary(gravmat);

    // compute synthetic gravity dsyn
    fvec dsyn;
    printf("forward computing ...\n");
    compute_gravity(vs,vsref,gmat,dsyn);
    if(remove_avg) {
        float mean = dsyn.sum() / dsyn.size();
        dsyn = dsyn - mean;
    }

    // write synthetic synthetic data
    printf("saving results ...\n");
    FILE *fp;
    infile.open(gravfile); 
    if(infile.fail()) {
        printf("cannot open %s\n",gravfile.data());
        exit(1);
    }
    fp = fopen(outfile.c_str(),"w");
    int np = dsyn.size();
    for(int i = 0; i < np; i ++){
        float a,b,g;
        infile >> a >> b >> g;
        fprintf(fp,"%g %g %g\n",a,b,dsyn[i]);
    }
    fclose(fp);
    infile.close();

    return 0;
}