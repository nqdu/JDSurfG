#pragma once
template<typename T>
class LSMRDict{
    public:
    T atol,btol;
    T conlim; 
    int itnlim,istop,itn;
    T anorm,acond,arnorm;
    T xnorm,rnorm;
    T damp,weight;
    int localSize;
    bool verbose;
    
    // some parameters defined
    // you could change it according to your own settings
    LSMRDict(int iterlim,int dimension,T damp0,T weight0){
        atol = 1.0e-5,btol = 1.0e-5;
        conlim = 1.0e6;
        istop = 0;

        anorm = 0.0,acond = 0.0,arnorm = 0.0;
        xnorm = 0.0,rnorm = 0.0;

        itnlim = iterlim;
        localSize = dimension;
        damp = damp0;
        weight = weight0;

        verbose = false;
    }
};


extern "C"{
void LSMR(int m, int n, int lenrw, int *rw, int *col, 
        float *val, float *b, float damp, float atol, float btol,
        float conlim, int itnlim, int localSize, float *x, 
        int *istop, int *itn, float *normA, float *condA, 
        float *normr, float *normAr, float *normx);

void LSMR_csr(int m, int n, float *val,int *indices,int *indptr, 
            float *b, float damp, float atol, float btol,
        float conlim, int itnlim, int localSize, float *x, 
        int *istop, int *itn, float *normA, float *condA, 
        float *normr, float *normAr, float *normx,int verbose);
}