#ifndef JDSURFG_SHARED_CSR_MATRIX_H_
#define JDSURFG_SHARED_CSR_MATRIX_H_

#include<iostream>

template<typename T=float>
class LSMRDict{
public:
    T atol,btol;
    T conlim; 
    int itnlim,istop,itn;
    T anorm,acond,arnorm;
    T xnorm,rnorm;
    T damp,weight;
    int localSize,num_threads;
    bool verbose;
    
    // some parameters defined
    // you could change it according to your own settings
    LSMRDict(int iterlim,int dimension,T damp0,T weight0,int nthreads){
        atol = 1.0e-5,btol = 1.0e-5;
        conlim = 1.0e6;
        istop = 0;

        anorm = 0.0,acond = 0.0,arnorm = 0.0;
        xnorm = 0.0,rnorm = 0.0;

        itnlim = iterlim;
        localSize = dimension;
        damp = damp0;
        weight = weight0;
        num_threads = nthreads;

        verbose = false;
    }
};


/**
 * Compressed Sparse Row matrix Class
 * Parameters:
 * ---------------------------------------
 *  for i-th row, the col number for nonzeros is indices[indptr[i]:indptr[i+1]]
 *  See Scipy doc for more information
 *  https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
 */
class csr_matrix{

private:
    int MATRIX_ROW,MATRIX_COL;
    bool ALLOCATE = false;

public:
    // row_index :indices
    int *indices,*indptr;
    float *data;
    int nonzeros;

    void initialize(int rows,int cols,int nar);

    // constructor
    csr_matrix(int rows,int cols,int nar);

    //contructor
    csr_matrix();

    // destructor
    ~csr_matrix();

    void aprod(int mode,float *x, float *y) const ;
    int rows() const;
    int cols() const;

    void cpp2fort(bool inverse=false);

    void lsmr_solver(float *x,const float *b,LSMRDict<float> &dict);

    void read_binary(const std::string &filename);
    void write_binary(const std::string &filename) const;
}; 

// merge binary files
void merge_csr_files(int nprocs,const std::string &outfile);
void add_regularization(csr_matrix &smat,float weight,int nx,int ny,int nz);
 

#endif // end JDSURFG_SHARED_CSR_MATRIX_H_
