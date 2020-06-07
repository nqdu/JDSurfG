#include<iostream>
#include<string>
#include<fstream>
#include"JointTomo.hpp"
using namespace Eigen;

void print_mean_and_rms(VectorXf a,std::string info)
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
        modtrue = modfile +".true";
        refmod = "None";
        if(argc == 7) refmod = argv[6];
        if(argc == 8){
            refmod = argv[6];
            modtrue = argv[7];
        }
    }
    else{
        std::cout <<"Please run this executable file by:"<<std::endl;
        std::string info = "./this paramfile surfdatafile ";
        info += "gravdatafile gravmat initmod (refmod) (truemod)"; 
        std::cout <<info << std::endl;
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
    VectorXf dv(tomo.n); // variation of S-wave velocity for every iteration
    VectorXf dsyn(tomo.surf.num_data); // synthetic traveltime for every iteration
    VectorXf dg(tomo.obsg.np); // synthetic gravity for every iteration
    Tensor<float,3> vsf = tomo.mod.vs * 1.0f; // set inversion model to initial 

   // checkerboard test if required
    if(tomo.param.ifsyn){
        std::cout<<" Checkerboard Resolution Test Begin ..." << std::endl;
        tomo.checkerboard();
    }

    // build a folder to store synthetics for every iteration
    int flag = system("mkdir -p results");

   // save initial model
    int nx = tomo.mod.nx, ny = tomo.mod.ny,nz=tomo.mod.nz;
    std::string resfile="results/mod_iter"+ std::to_string(0) +  ".dat";
    outfile.open(resfile);
    for(int k=0;k<nz-1;k++){
    for(int j=0;j<ny-2;j++){
    for(int i=0;i<nx-2;i++){
        outfile << tomo.mod.lon(j+1) << " " <<  tomo.mod.lat(i+1) << " " 
                << tomo.mod.dep(k) << " " << vsf(i+1,j+1,k) << std::endl;
    }}}
    outfile.close();

    // inversion begin
    Map<VectorXf> Gr(tomo.obsg.Gr,tomo.obsg.np);
    int nt = tomo.surf.num_data;
    int ng = tomo.obsg.np;
    std::string info;
    for(int iter =0;iter < tomo.param.maxiter;iter++){
        std::cout << std::endl;
        std::cout << "Iteration " + std::to_string(iter+1) << std::endl;
        tomo.inversion(vsf,dv,dsyn,dg);

        // compute mean and rms of residuals
        VectorXf res(tomo.num_data);
        res.segment(0,nt) = tomo.surf.obst - dsyn;
        res.segment(nt,ng) = Gr - dg;

        // for traveltime
        info = "mean and rms of traveltime resuduals before this iteration(s):";
        print_mean_and_rms(res.segment(0,nt),info);
        
        // for gravity
        info="mean and rms of gravity resuduals before this iteration(mGal):";
        print_mean_and_rms(res.segment(nt,ng),info);

        // for joint residuals
        info="mean and rms of joint resuduals before this iteration:";
        print_mean_and_rms(res,info);

        // save current synthetics
        std::string resfile = "results/res"+std ::to_string(iter)+".dat";
        outfile.open(resfile);
        for(int i=0;i<tomo.num_data;i++){
            outfile << tomo.surf.dist(i) << " " << tomo.surf.obst(i) << " "\
                   << dsyn(i) <<std::endl;
        }
        outfile.close();

        resfile = "results/res_grav"+std ::to_string(iter)+".dat";
        outfile.open(resfile);
        for(int i=0;i<tomo.obsg.np;i++){
            outfile << Gr(i) << " "<< dg(i) <<std::endl;
        }
        outfile.close();

        // save current model
        resfile="results/mod_iter"+std ::to_string(iter+1)+".dat";
        outfile.open(resfile);
        for(int k=0;k<nz-1;k++){
        for(int j=0;j<ny-2;j++){
        for(int i=0;i<nx-2;i++){
            outfile << tomo.mod.lon(j+1) << " " <<  tomo.mod.lat(i+1) << " " 
                    << tomo.mod.dep(k) << " " << vsf(i+1,j+1,k) << std::endl;
        }}}
        outfile.close();
    }

    // compute theoretical traveltimes for last iteration
    std::cout <<std::endl;
    std::cout << " synthetic traveltime for the result model " << std::endl;
    tomo.forward(vsf,dsyn,dg);

    // compute mean and rms of residuals
    VectorXf res(tomo.num_data);
    res.segment(0,nt) = tomo.surf.obst - dsyn;
    res.segment(nt,ng) = Gr - dg;
    // for traveltime
    info = "mean and rms of traveltime resuduals after inversion(s):";
    print_mean_and_rms(res.segment(0,nt),info);
    
    // for gravity
    info="mean and rms of gravity resuduals after inversion (mGal):";
    print_mean_and_rms(res.segment(nt,ng),info);

    // for joint residuals
    info="mean and rms of joint resuduals after inversion:";
    print_mean_and_rms(res,info);
    std::cout << std::endl;

    std::cout << "Program finishes Successfully!" << std::endl;
    return 0;
}