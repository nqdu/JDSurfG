#include"tomography.hpp"
#include<fstream>
using Eigen::Tensor;
using Eigen::VectorXf;

void print_mean_and_rms(VectorXf a,std::string info)
{
    float mean = a.mean();
    float rms = a.array().pow(2.).mean() - mean * mean;
    rms = sqrt(rms);
    std::cout << info << " " << mean << " " << rms << std::endl;
}

int main(int argc, char* argv[]){
   // check input parameters
    std::string paramfile,modfile,datafile,modtrue;
    std::cout << "\n";
    std::cout << "\t Direct Surface Wave Tomography" << std::endl;
    std::cout<<std::endl;
    if(argc == 1 ){
        std::cout <<"NO INPUT FILES ARE given!, now use the default ones ..."<<std::endl;
        paramfile="ex/DSurfTomo.in"; 
        modfile="ex/MOD";
        datafile="ex/surfdataTB.dat";
        modtrue="ex/MOD.true";
        std::cout <<"default files are: "<<paramfile << " " << datafile
                  << " " << modfile <<" " << modtrue <<std::endl;
    }
    else if(argc == 4 || argc == 5){
        paramfile = argv[1];
        modfile = argv[3];
        datafile = argv[2];
        modtrue = "None";
        if(argc == 5) modtrue = argv[4];
    }
    else if(argc == 2 && !strcmp(argv[1],"-h")){
        std::cout <<"Please run this executable file by:"<<std::endl;
        std::cout <<"./this paramfile datafile initmod (truemod)" << std::endl;
        exit(0);
    }
    else{
        std::cout <<"Please run this executable file -h for help"<<std::endl;
        exit(0);
    }
    std::cout <<std::endl;

    // output stream
    std::ofstream outfile;

    // class for inversion;
    DSurfTomo tomo;

    // read all the parameters, initial model, and observed data
    tomo.readdata(paramfile,datafile,modfile,modtrue);

   // initialize some parameters
    VectorXf dsyn(tomo.num_data); // synthetics for every step
    Tensor<float,3> vsf = tomo.mod.vs * 1.0f;  // set inversion model to initial one

   // checkerboard test if required
    if(tomo.param.ifsyn == 1){
        std::cout<<" Checkerboard Resolution Test Begin ..." << std::endl;
        tomo.checkerboard();
    }

    // build a folder to store synthetics for every iteration
    int flag = system("mkdir -p results");

    // save initial model
    int nx = tomo.mod.nx, ny = tomo.mod.ny,nz=tomo.mod.nz;
    std::string resfile="results/mod_iter"+ std::to_string(0) +  ".dat";
    outfile.open(resfile);
    for(int k=0;k<nz;k++){
    for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++){
        outfile << tomo.mod.lon(j) << " " <<  tomo.mod.lat(i) << " " 
                << tomo.mod.dep(k) << " " << vsf(i,j,k) << std::endl;
    }}}
    outfile.close();

    // inversion begin
    std::string info;
    for(int iter =0;iter < tomo.param.maxiter;iter++){
        std::cout << std::endl;
        std::cout << "Iteration " + std::to_string(iter+1) << std::endl;
        tomo.inversion(vsf,dsyn);

        // compute mean and rms of residuals
        VectorXf res = tomo.surf.obst - dsyn;
        info = "mean and rms before this iteration(s):";
        print_mean_and_rms(res,info);

        // save current synthetics
        std::string resfile = "results/res"+std ::to_string(iter)+".dat";
        tomo.surf.write_disper(dsyn,resfile);

        // save current model
        resfile="results/mod_iter"+std ::to_string(iter+1)+".dat";
        outfile.open(resfile);
        for(int k=0;k<nz;k++){
        for(int j=0;j<ny;j++){
        for(int i=0;i<nx;i++){
            outfile << tomo.mod.lon(j) << " " <<  tomo.mod.lat(i) << " " 
                    << tomo.mod.dep(k) << " " << vsf(i,j,k) << std::endl;
        }}}
        outfile.close();
    }

    // compute theoretical traveltimes for last iteration
    std::cout <<std::endl;
    std::cout << " synthetic traveltime for the result model " << std::endl;
    tomo.forward(vsf,dsyn);

    // compute mean and rms of residuals
    VectorXf res = tomo.surf.obst - dsyn;
    info = "mean and rms of traveltime residuals before this iteration(s):";
    print_mean_and_rms(res,info);

    // save synthetics for last iteration
    int maxiter = tomo.param.maxiter;
    resfile = "results/res"+std ::to_string(maxiter)+".dat";
    tomo.surf.write_disper(dsyn,resfile);

    std::cout << std::endl;
    int ierr = system("rm -r kernel");
    std::cout << "Program finishes Successfully!" << std::endl;
    return 0;

}