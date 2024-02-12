#include "invparam.hpp"
#include "IOFunc.hpp"
#include <sstream>

void InverseParamsBase :: 
read_file(const std::string &paramfile) {
    // open file
    std::ifstream infile; infile.open(paramfile);

    // inv method
    read_par_regex("INV_METHOD",inv_method,infile);

    // read parameters from file
    read_par_regex("NITERS",maxiter,infile);
    int ierr = read_par_regex("ITER_CURRENT",iter_cur,infile);
    if(ierr == 1) {
        iter_cur = 0;
    }

    // constraints
    read_par_regex("MIN_VELOC",minvel,infile);
    read_par_regex("MAX_VELOC",maxvel,infile);

    // read inv params based on inv_method
    if(inv_method == 0) { //LSMR
        read_par_regex("SMOOTH",smooth,infile);
        read_par_regex("DAMP",damp,infile);
        read_par_regex("NTHREADS",nthreads,infile);
    } 
    else {
        read_par_regex("SMOOTH_IN_KM",smooth_in_km,infile);
        read_par_regex("SIGMA_H",sigma_h,infile);
        read_par_regex("SIGMA_V",sigma_v,infile);
        ierr = read_par_regex("ITER_START",iter_start,infile);
        if(ierr == 1) {
            iter_start = 0;
        }
    }

    // synthetic test 
    read_par_regex("SYN_TEST",ifsyn,infile);

    // print on the screen 
    printf("Inversion Parameters:\n");
    printf("===================================\n");
    printf("Min amd max velocity(km/s) = %f, %f\n",minvel,maxvel);
    printf("Max iterations = %d\n",maxiter);
    printf("current model = %d\n",iter_cur);

    if(inv_method == 0) {
        printf("use LSMR solver: ");
        printf("Number of Threads Used = %d\n",nthreads);
        printf("smooth = %f, damp = %f\n",smooth,damp);
    }
    else {
        if(inv_method == 1) {
            printf("use nonlinear-CG optimization:\n");
        }
        else {
            printf("use LBFGS optimization:\n");
        }
        printf("sigma_h = %f,  sigma_v = %f\n",sigma_h,sigma_v);

        // read line search params
        read_par_regex("MAX_REL_STEP",MAX_REL_STEP,infile);
    }

    // close file
    infile.close();
    
}