#include "surftomo/surftomo.hpp"
#include <fstream>
#include "shared/IOFunc.hpp"
#include "nonlinear_opt.hpp"

int main(int argc, char* argv[]){
   // check input parameters
    std::string paramfile,modfile,datafile,modtrue;
    printf("\n**************************************\n");
    printf("*** Direct Surface Wave Tomography ***\n");
    printf("**************************************\n");
    if(argc == 4 || argc == 5){
        paramfile = argv[1];
        modfile = argv[3];
        datafile = argv[2];
        modtrue = "None";
        if(argc == 5) modtrue = argv[4];
    }
    else {
        printf("Please run this executable file by:\n");
        printf("./this paramfile datafile initmod (truemod)\n");
        exit(1);
    }
    printf("\n");

    // class for inversion;
    DSurfTomo tomo;

    // read all the parameters, initial model, and observed data
    tomo.read_invparams(paramfile);
    tomo.read_model(modfile,modtrue);
    tomo.read_data(datafile);
    const auto &param = tomo.param;
    
   // initialize some parameters
    fvec dsyn(tomo.surf.obst.size()); // synthetics for every step
    fmat3 vsf = tomo.vsinit * 1.0f;  // set inversion model to initial one

   // checkerboard test if required
    if(tomo.param.ifsyn == 1){
        printf("\nCheckerboard Resolution Test Begin ...\n" );
        tomo.checkerboard();
    }

    // build a folder to store synthetics for every iteration
    std::string outdir = "results";
    create_directory(outdir.c_str());

    // save initial model
    int nx = vsf.dimension(0), ny = vsf.dimension(1);
    int nz= vsf.dimension(2);
    int n = (nx-2) * (ny-2) * (nz-1);
    std::string resfile= outdir + "/mod_iter"+ std::to_string(param.iter_cur) +  ".dat";
    FILE *fp = fopen(resfile.c_str(),"w");
    for(int k=0;k<nz;k++){
    for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++){
        fprintf(fp,"%g %g %g %g\n",tomo.lon[j],tomo.lat[i],tomo.dep[k],vsf(i,j,k));
    }}}
    fclose(fp);

    // allocate space for input x 
    fvec x;
    if(tomo.param.inv_method > 0){
        x.resize(n);
        int ic = 0;
        for(int k=0;k<nz-1;k++){
        for(int j=0;j<ny-2;j++){
        for(int i=0;i<nx-2;i++){
            x[ic] = vsf(i+1,j+1,k);
            ic += 1;
        }}}
    }

    // inversion begin
    std::string info;
    NonlinOPT<DSurfTomo> opt(outdir,param.MAX_REL_STEP,param.iter_start,param.iter_cur);
    for(int ii = 0; ii < param.maxiter;ii ++){
        int iter = ii + param.iter_cur;
        printf("\n");
        printf("Iteration %d\n",iter+1);
        if(param.inv_method == 0) {
            tomo.inversion(vsf,dsyn);
        }
        else {
            std::string method = "CG";
            if(param.inv_method == 2) {
                method = "LBFGS";
            }
            opt.update(tomo,x,dsyn,method);

            // update vsf
            int ic = 0;
            for(int k=0;k<nz-1;k++){
            for(int j=0;j<ny-2;j++){
            for(int i=0;i<nx-2;i++){
                vsf(i+1,j+1,k) = x[ic];
                ic += 1;
            }}}
        }

        // compute mean and rms of residuals
        fvec res = tomo.surf.obst - dsyn;
        float mean = res.sum() / tomo.surf.obst.size();
        float rms = tomo.compute_misfit(dsyn) / res.size();
        printf("mean and rms before this iteration(s): %g %g\n",mean,rms);

        // save current synthetics
        std::string resfile = outdir + "/res"+std ::to_string(iter)+".dat";
        tomo.surf.write_syn(dsyn,resfile);

        // save current model
        resfile = outdir + "/mod_iter"+std ::to_string(iter+1)+".dat";
        fp = fopen(resfile.c_str(),"w");
        for(int k=0;k<nz;k++){
        for(int j=0;j<ny;j++){
        for(int i=0;i<nx;i++){
            fprintf(fp,"%g %g %g %g\n",tomo.lon[j],tomo.lat[i],tomo.dep[k],vsf(i,j,k));
        }}}
        fclose(fp);
    }

    // compute theoretical traveltimes for last iteration
    printf("\n");
    printf("synthetic traveltime for the result model \n");
    tomo.surf.travel_time(vsf,dsyn);

    // compute mean and rms of residuals
    fvec res = tomo.surf.obst - dsyn;
    float mean = res.sum() / tomo.surf.obst.size();
    float rms = tomo.compute_misfit(dsyn) / res.size();
    printf("mean and rms before this iteration(s): %g %g\n",mean,rms);

    // save synthetics for last iteration
    int maxiter = tomo.param.maxiter + param.iter_cur;;
    resfile = outdir + "/res"+std ::to_string(maxiter)+".dat";
    tomo.surf.write_syn(dsyn,resfile);

    printf("\n");
    printf("Program finishes Successfully!\n\n");
    
    return 0;

}