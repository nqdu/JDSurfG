#include"gravmat.hpp"
#include"defmod.hpp"
#include"coo_matrix.hpp"
#include<fstream>

int main(int argc,char *argv[]){
    // check input parameters
    std::string paramfile,modfile,datafile,modtrue;
    if(argc == 1 ){
        std::cout <<"NO INPUT FILES ARE given!, now use the default ones ..."<<std::endl;
        paramfile="DSurfTomo.in"; 
        datafile="obsgrav.dat";
        modfile="MOD";
        modtrue="MOD.true";
        std::cout <<"default files are: "<<paramfile << " " << datafile
                  << " " << modfile <<" " << modtrue <<std::endl;
    }
    else if(argc == 4 || argc == 5){
        paramfile = argv[1];
        modfile = argv[3];
        datafile = argv[2];
        modtrue = modfile +".true";
        if(argc == 5) modtrue = argv[4];
    }
    else{
        std::cout <<"Please run this executable file by:"<<std::endl;
        std::cout <<"./this paramfile datafile initmod (truemod)" << std::endl;
        exit(0);
    }
    std::cout <<std::endl;
    MOD3DSphGra mod;
    OBSSphGraRandom obs;

    // read density model and gravity data
    mod.read_model(paramfile,modfile,modtrue);
    obs.read_obs_data(datafile);

    // init sparse matrix
    int n = mod.nx * mod.ny * mod.nz;
    int nonzeros = (int)(n * obs.np * 0.1);
    coo_matrix<float> smat(obs.np,n,nonzeros);

    //change coordinates
    obs.chancoor(1);
    mod.chancoor(1);

    // compute gravity matrix and synthetic gravity data
    FILE *fp;
    if(mod.synflag){
        std::cout << "begin to compute gravity matrix and also synthesize gravity data"
                  << std::endl;
        gravmat_parallel(mod,obs,smat);
        float *x = new float [n];
        float *y = new float [obs.np];
        for(int i=0;i<n;i++){
            x[i] = mod.density[i] - mod.density0[i];
        }
        smat.aprod(1,x,y); 

        // save gravity data
        fp = fopen("gravity.dat","w");
        if(fp == NULL){
            std::cout <<"cannot open file gravity.dat";
            exit(0);
        }
        obs.chancoor(0);
        for(int i=0;i<obs.np;i++){
            int flag = fprintf(fp,"%f %f %g\n",obs.lon[i],obs.lat[i],y[i]);
            if(flag <0){
                printf("cannot write in gravity.dat\n");
                exit(0);
            }
        }
        fclose(fp);

        // free space
        delete[] x;
        delete[] y;
    }
    else{
        std::cout << "begin to compute gravity matrix" << std::endl;
        gravmat_parallel(mod,obs,smat);
    }

    // save gravity matrix
    fp = fopen("gravmat.dat","w");
    if(fp == NULL){
        std::cout <<"cannot open file gravmat.dat";
        exit(0);
    }
    for(int i=0;i<smat.nonzeros;i++){
        fprintf(fp,"%d %d %g\n",smat.rw[i],smat.col[i],smat.val[i]);
    }
    
    fclose(fp);

}