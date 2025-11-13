#include "invparam.hpp"
#include "IOFunc.hpp"
#include <sstream>

/**
 * @brief read inversion parameters from file
 * 
 * @param paramfile parameter file name
 */
void InverseParamsBase :: 
read_file(const std::string &paramfile) {
    // open file
    std::ifstream infile; infile.open(paramfile);

    // macros to read inversion parameters
    #define READ_PARAM(NAME,VAR) \
        ierr = read_par_regex(NAME,VAR,infile); \
        if(ierr == 1){ \
            printf("cannot find %s\n",NAME); \
            exit(1); \
        }

    // inv method
    int ierr;
    READ_PARAM("INV_METHOD",inv_method);

    // read parameters from file
    READ_PARAM("NITERS",maxiter);
    ierr = read_par_regex("ITER_CURRENT",iter_cur,infile);
    if(ierr == 1) {
        iter_cur = 0;
    }

    // constraints
    READ_PARAM("MIN_VELOC",minvel);
    READ_PARAM("MAX_VELOC",maxvel);

    // read inv params based on inv_method
    if(inv_method == 0) { //LSMR
        READ_PARAM("SMOOTH",smooth);
        READ_PARAM("DAMP",damp);
        READ_PARAM("NTHREADS",nthreads);
    } 
    else {
        READ_PARAM("SMOOTH_IN_KM",smooth_in_km);
        READ_PARAM("SIGMA_H",sigma_h);
        READ_PARAM("SIGMA_V",sigma_v);
        ierr = read_par_regex("ITER_START",iter_start,infile);
        if(ierr == 1) {
            iter_start = 0;
        }
    }

    // synthetic test 
    READ_PARAM("SYN_TEST",ifsyn);

    // read topography if required
    ierr = read_par_regex("TOPO_CORR",topo_corr,infile);
    if(ierr == 1) {
        topo_corr = 0;
    }

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
        READ_PARAM("MAX_REL_STEP",MAX_REL_STEP);
    }

    // close file
    infile.close();

    #undef READ_PARAM
    
}