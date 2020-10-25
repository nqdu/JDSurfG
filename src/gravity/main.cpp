#include"gravmat.hpp"
#include<fstream>
#include<string.h>

int main(int argc,char *argv[]){
    // check input parameters
    std::string paramfile,datafile,modfile;
    int ierr;
    if(argc == 1 ){
        std::cout <<"NO INPUT FILES ARE given!, now use the default ones ..."<<std::endl;
        paramfile="DSurfTomo.in"; 
        datafile="obsgrav.dat";
        modfile="MOD";
        std::cout <<"default files are: "<<paramfile << " " << datafile
                  << " " << modfile <<std::endl;
    }
    else if(argc == 4){
        paramfile = argv[1];
        modfile = argv[3];
        datafile = argv[2];
    }
    else if (argc == 2 && !strcmp(argv[1],"-h")){
        std::cout <<"Please run this executable file by:"<<std::endl;
        std::cout <<"./this paramfile datafile vs-model " << std::endl;
        exit(0);
    }
    else{
        std::cout <<"Please run this executable file -h:"<<std::endl;
        exit(0);
    }
    std::cout <<std::endl;
    MOD3DSphGra mod;
    OBSSphGraRandom obs;

    // read density model and gravity data
    //mod.read_model(paramfile,modfile);
    mod.read_model(paramfile,modfile);
    obs.read_obs_data(datafile);

    // init sparse matrix
    int n = mod.nx * mod.ny * mod.nz;
    int nonzeros = (int)(0.1 * n * obs.np);
    csr_matrix<float> smat(obs.np,n,nonzeros);

    //change coordinates
    obs.chancoor(1);
    mod.chancoor(1);

    // compute gravity matrix and synthetic gravity data
    FILE *fp;
    std::cout << "begin to compute gravity matrix" << std::endl;
    gravmat_parallel(mod,obs,smat); 

    // save gravity matrix
    fp = fopen("gravmat.dat","w");
    if(fp == NULL){
        std::cout <<"cannot open file gravmat.dat";
        exit(0);
    }
    for(int i=0;i<smat.rows();i++){
        int start = smat.indptr[i];
        int end = smat.indptr[i+1];
        int nar = end - start;
        ierr = fprintf(fp,"# %d %d\n",i,nar);
        for(int j=start;j<end;j++){
            ierr = fprintf(fp,"%d %g\n",smat.indices[j],smat.data[j]);
        }
    }
    
    fclose(fp);

}