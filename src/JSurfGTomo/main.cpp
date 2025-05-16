#include "JSurfGTomo/JSurfGTomo.hpp"
#include <fstream>
#include "shared/IOFunc.hpp"
#include "nonlinear_opt.hpp"

int main(int argc, char* argv[]){
   // check input parameters
    std::string paramfile,modfile,swdfile,gravfile,gravmat,modtrue,modref;
    printf("\n**************************************\n");
    printf("*** Joint Inversion of SWD and Gravity ***\n");
    printf("**************************************\n");
    if(argc == 6 || argc == 7 || argc == 8){
        paramfile = argv[1];
        swdfile = argv[2];
        gravfile = argv[3];
        gravmat = argv[4];
        modfile = argv[5];
        modtrue = "None";
        modref = "None";
        if(argc == 7) modref = argv[6];
        if(argc == 8) {
            modref = argv[6];
            modtrue = argv[7];
        }
    }
    else {
        printf("Please run this executable file by:\n");
        printf("./this paramfile swdfile gravfile gravmat initmod (refmod) (truemod)\n");
        exit(1);
    }
    printf("\n");

    // class for inversion;
    JSurfGTomo tomo;

    // read all the parameters, initial model, and observed data
    tomo.read_invparams(paramfile);
    tomo.read_model(modfile,modtrue,modref);
    tomo.read_data(swdfile,gravfile);
    tomo.read_gravmat(gravmat);
    const auto &param = tomo.param;

    // data dimen
    int m1 = tomo.surf.obst.size(), m2 = tomo.obsg.size();
    fvec res1(m1),res2(m2);
    
   // initialize some parameters
    fvec dsyn(m1+m2); // synthetics for every step
    fmat3 vsf = tomo.vsinit * 1.0f;  // set inversion model to initial one

   // checkerboard test if required
    if(tomo.param.ifsyn == 1){
        printf("\nCheckerboard Resolution Test Begin ...\n" );
        tomo.checkerboard();
    }

    // build a folder to store synthetics for every iteration
    std::string outdir = "resultsJ";
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
    NonlinOPT <JSurfGTomo> opt(outdir,param.MAX_REL_STEP,param.iter_start,param.iter_cur);
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
        res1 = tomo.surf.obst - dsyn.segment(0,m1);
        res2 = tomo.obsg - dsyn.segment(m1,m2);
        float mean1 = res1.sum() / m1, mean2 = res2.sum() / m2;
        float rms1 = std::sqrt(res1.square().sum() / m1);
        float rms2 = std::sqrt(res2.square().sum() / m2);;
        printf("mean and rms of SWD before this iteration(s): %g %g\n",mean1,rms1);
        printf("mean and rms of Gravity before this iteration(s): %g %g\n",mean2,rms2);
        printf("joint misfit = %g\n",tomo.compute_misfit(dsyn));

        // save current synthetics
        std::string resfile1 = outdir + "/res_swd"+std ::to_string(iter)+".dat";
        std::string resfile2 = outdir + "/res_grav"+std ::to_string(iter)+".dat";
        tomo.write_syn(dsyn,resfile1,resfile2);

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
    fvec dsyn1(m1),dsyn2(m2);
    tomo.surf.travel_time(vsf,dsyn1);
    tomo.compute_gravity(vsf,dsyn2);

    // compute mean and rms of residuals
    res1 = tomo.surf.obst - dsyn1;
    res2 = tomo.obsg - dsyn2;
    float mean1 = res1.sum() / m1, mean2 = res2.sum() / m2;
    float rms1 = std::sqrt(res1.square().sum() / m1);
    float rms2 = std::sqrt(res2.square().sum() / m2);
    printf("mean and rms of SWD of final model: %g %g\n",mean1,rms1);
    printf("mean and rms of Gravity of final model: %g %g\n",mean2,rms2);
    printf("joint misfit = %g\n",tomo.compute_misfit(dsyn));

    // save synthetics for last iteration
    int maxiter = tomo.param.maxiter + param.iter_cur;
    std::string resfile1 = outdir + "/res_swd"+std ::to_string(maxiter)+".dat";
    std::string resfile2 = outdir + "/res_grav"+std ::to_string(maxiter)+".dat";
    tomo.write_syn(dsyn,resfile1,resfile2);

    printf("\n");
    printf("Program finishes Successfully!\n\n");
    
    return 0;

}