#include "gravity_module.hpp"
#include<fstream>
#include<string.h>
#include "shared/csr_matrix.hpp"

int main(int argc,char *argv[]){

    // check input parameters
    std::string modfile,datafile;
    if(argc == 3){
        modfile = argv[1];
        datafile = argv[2];
    }
    else {
        printf("Please run this executable file by:\n");
        printf("./mkmat vs-model datafile  \n");
        printf("Example: ./mkmat MOD obsgrav.dat\n");
        exit(1);
    }

    printf("\n");
    MOD3DSphGra mod;
    OBSSphGraRandom obs;

    // read density model and gravity data
    //mod.read_model(paramfile,modfile);
    mod.read_model(modfile);
    obs.read_obs_data(datafile);

    //change coordinates
    obs.chancoor(1);
    mod.chancoor(1);

    // compute gravity matrix and synthetic gravity data
    printf("\nbegin to compute gravity matrix ...\n");
    generate_gravmat(mod,obs,"gravmat.bin");
    printf("\n");

    return 0;
}