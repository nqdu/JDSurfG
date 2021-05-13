#include<iostream>
#include<string>
#include<fstream>
#include"tomography.hpp"
using Eigen::VectorXf;
using Eigen::Tensor;

void print_mean_and_rms(const VectorXf &a,std::string info,float weight1=0.0,float weight2=0.0)
{
    float mean = a.mean();
    float rms = a.array().pow(2.).mean() - mean * mean;
    rms = sqrt(rms);
    std::cout << info << " " << mean << " " << rms << std::endl;
}

int main(int argc, char* argv[]){
    // check input parameters
    std::string paramfile,modfile,surfdata,modtrue;
    std::string gravmat,gravdata,refmod;
    std::cout << "\n";
    std::cout << "\t Joint Inversion of Direct Surface Wave Tomography and Gravity" << std::endl;
    std::cout<<std::endl;
    if(argc == 1 ){
        std::cout <<"NO INPUT FILES ARE given!, now use the default ones ..."<<std::endl;
        paramfile="JointSG.in"; 
        modfile="MOD.surf";
        surfdata="surfdataSC.new.dat";
        modtrue="MOD.true";
        gravmat="gravmat.dat";
        gravdata="gravity_obs.dat";
        refmod="MOD";
        std::cout <<"default files are: "<<paramfile <<
                  " " << surfdata << " " << gravdata <<
                  " " << gravmat << " " << modfile << " " 
                  << refmod << " " << modtrue<<std::endl;
    }
    else if(argc == 6 || argc == 7 || argc == 8){
        paramfile = argv[1];
        surfdata = argv[2];
        gravdata = argv[3];
        gravmat = argv[4];
        modfile = argv[5];
        modtrue = "None";
        refmod = "None";
        if(argc == 7) refmod = argv[6];
        if(argc == 8){
            refmod = argv[6];
            modtrue = argv[7];
        }
    }
    else if(argc == 2 && !strcmp(argv[1],"-h")){
        std::cout <<"Please run this executable file by:"<<std::endl;
        std::string info = "./this paramfile surfdatafile ";
        info += "gravdatafile gravmat initmod (refmod) (truemod)"; 
        std::cout <<info << std::endl;
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
    JointTomo tomo;

    // read all the parameters, initial model, and observed data
    tomo.readdata(paramfile,modfile,surfdata,gravdata,gravmat,refmod,modtrue);

    // initialize some parameters
    VectorXf dsyn(tomo.surf.num_data); // synthetic traveltime for every iteration
    VectorXf dg(tomo.obsg.np); // synthetic gravity for every iteration
    Tensor<float,3> vsf = tomo.mod.vs * 1.0f; // set inversion model to initial 

   // checkerboard test if required
    if(tomo.param.ifsyn==1){
        std::cout<<"Checkerboard Resolution Test Begin ..." << std::endl;
        tomo.checkerboard();
    }

    // build a folder to store synthetics for every iteration
    int flag = system("mkdir -p results");

   // save initial model
    int nx = tomo.mod.nx, ny = tomo.mod.ny,nz=tomo.mod.nz;
    std::string resfile="results/joint_mod_iter"+ std::to_string(0) +  ".dat";
    outfile.open(resfile);
    for(int k=0;k<nz;k++){
    for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++){
        outfile << tomo.mod.lon(j) << " " <<  tomo.mod.lat(i) << " " 
                << tomo.mod.dep(k) << " " << vsf(i,j,k) << std::endl;
    }}}
    outfile.close();

    // remove mean of gravity data
    Eigen::Map<VectorXf> Gr(tomo.obsg.Gr,tomo.obsg.np);
    float mean = Gr.sum() / tomo.obsg.np;
    Gr.array() -= mean;
    
    // inversion begin
    int nt = tomo.surf.num_data;
    int ng = tomo.obsg.np;
    std::string info;
    for(int iter =0;iter < tomo.param.maxiter;iter++){
        std::cout << std::endl;
        std::cout << "Iteration " + std::to_string(iter+1) << std::endl;
        tomo.inversion(vsf,dsyn,dg);

        // compute mean and rms of residuals
        VectorXf res(tomo.num_data);
        res.segment(0,nt) = tomo.surf.obst - dsyn;
        res.segment(nt,ng) = Gr - dg;

        // for traveltime
        info = "mean and rms of traveltime residuals before this iteration(s):";
        print_mean_and_rms(res.segment(0,nt),info);
        
        // for gravity
        info="mean and rms of gravity residuals before this iteration(mGal):";
        print_mean_and_rms(res.segment(nt,ng),info);

        // for joint residuals
        info="mean and rms of joint residuals before this iteration:";
        print_mean_and_rms(res,info);

        // save current synthetics
        std::string resfile = "results/res_surf"+std ::to_string(iter)+".dat";
        tomo.surf.write_disper(dsyn,resfile);

        resfile = "results/res_grav"+std ::to_string(iter)+".dat";
        outfile.open(resfile);
        for(int i=0;i<ng;i++){
            outfile << tomo.obsg.lon[i] << " " << tomo.obsg.lat[i] << " " <<Gr(i) << " "<< dg(i) <<std::endl;
        }
        outfile.close();

        // save current model
        resfile="results/joint_mod_iter"+std ::to_string(iter+1)+".dat";
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
    tomo.forward(vsf,dsyn,dg);
    mean = dg.sum() / dg.size();
    dg.array() -= mean;

    // compute mean and rms of residuals
    VectorXf res(tomo.num_data);
    res.segment(0,nt) = tomo.surf.obst - dsyn;
    res.segment(nt,ng) = Gr - dg;
    // for traveltime
    info = "mean and rms of traveltime residuals after inversion(s):";
    print_mean_and_rms(res.segment(0,nt),info);
    
    // for gravity
    info="mean and rms of gravity residuals after inversion (mGal):";
    print_mean_and_rms(res.segment(nt,ng),info);

    // for joint residuals
    info="mean and rms of joint residuals after inversion:";
    print_mean_and_rms(res,info);
    std::cout << std::endl;

    // save synthetics for last iteration
    int maxiter = tomo.param.maxiter;
    resfile = "results/res_surf"+std ::to_string(maxiter)+".dat";
    tomo.surf.write_disper(dsyn,resfile);

    resfile = "results/res_grav"+std ::to_string(maxiter)+".dat";
    outfile.open(resfile);
    for(int i=0;i<ng;i++){
        outfile << tomo.obsg.lon[i] << " " << tomo.obsg.lat[i] << " " <<    Gr(i) << " "<< dg(i) <<std::endl;
    }
    outfile.close();

    int ierr = system("rm -r kernelJ");
    std::cout << "Program finishes Successfully!" << std::endl;

    return 0;
}
