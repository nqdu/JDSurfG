#pragma once
#include<iostream>
#include"lsmr.hpp"

/**
 * Compressed Sparse Row matrix Class
 * Parameters:
 * ---------------------------------------
 *  for i-th row, the col number for nonzeros is indices[indptr[i]:indptr[i+1]]
 *  See Scipy doc for more information
 *  https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
 */
template<typename T>
class csr_matrix{
    private:
    unsigned int MATRIX_ROW,MATRIX_COL;
    bool ALLOCATE;

    public:
    // row_index :indices
    int * __restrict__ indices,* __restrict__ indptr;
    T * __restrict__ data;
    unsigned int nonzeros;

    void initialize(int rows,int cols,int nar)
    {
        MATRIX_ROW = rows;
        MATRIX_COL = cols;
        nonzeros = nar;
        indices = new int [nar]();
        data = new T [nar]();
        indptr = new int[rows + 1]();
        ALLOCATE = true;
    }

    // constructor
    csr_matrix(unsigned int rows,unsigned int cols,unsigned int nar)
    {
        initialize(rows,cols,nar);
    }

    //contructor
    csr_matrix()
    {
        ALLOCATE = false;
    }

    // destructor
    ~csr_matrix(){
        if(ALLOCATE){
            delete [] data;
            delete [] indptr;
            delete[] indices;
            ALLOCATE = false;
        }
    }

    // get rows
    unsigned int rows(){
        return MATRIX_ROW;
    }

    unsigned int cols(){
        return MATRIX_COL;
    }

    bool is_alloc(){
        return ALLOCATE;
    }

    int get_nonzeros(){
        return nonzeros;
    }

    /* 
        matrix-vector multiplication
    Parameters:
    ------------------------------------------
    mode : 1.   compute y = y + A * x
        2.   compute x = x + A.T * y
    */
    void aprod(int mode,T *x, T *y){
        if(mode == 1){ // y = y + A * x
            for(int i=0;i<MATRIX_ROW;i++){
            for(int j=indptr[i];j<indptr[i+1];j++){
                y[i] += data[j] * x[indices[j]];
            }}
        }
        else{
            for(int i=0;i<MATRIX_ROW;i++){
            for(int j=indptr[i];j<indptr[i+1];j++){
                x[indices[j]] += data[j] * y[i];
            }}
        }
    }

    void cpp2fortran(bool inverse=false){
        if(!inverse){
            for(int i=0;i<nonzeros;i++)indices[i] +=1;
            for(int i=0;i<MATRIX_ROW+1;i++) indptr[i] +=1;
        }
        else{
            for(int i=0;i<nonzeros;i++)indices[i] -=1;
            for(int i=0;i<MATRIX_ROW+1;i++) indptr[i] -=1;
        }
    }

    void LsmrSolver(T *b,T *x,LSMRDict<float> &dict){

        // initialize x
        for(int i=0;i<MATRIX_COL;i++) x[i] = 0.0;
        cpp2fortran();

        // solve by lsmr module
        int show = 1;
        if(dict.verbose == false) show = 0;

        LSMR_csr(MATRIX_ROW,MATRIX_COL,data,indices,indptr,b,dict.damp,dict.atol,dict.btol,
            dict.conlim,dict.itnlim,dict.localSize,x,&dict.istop,
            &dict.itn,&dict.anorm,&dict.acond,&dict.rnorm,
            &dict.arnorm,&dict.xnorm,show);
        cpp2fortran(true);
    }

    void read(std::string filename){
        // open file
        FILE *fp;
        if((fp=fopen(filename.c_str(),"r"))==NULL){
            std::cout << "cannot open file "<< filename << std::endl;
            exit(0);
        } 

        // read from file
        char line[100],dummy;
        int rw_idx,nar;
        int ierr;
        
        while(fgets(line,sizeof(line),fp)!=NULL){
            if(line[0] != '#') break;
            sscanf(line,"%c%d%d",&dummy,&rw_idx,&nar);
            int start = indptr[rw_idx];
            indptr[rw_idx + 1] = nar + start;
            int end = indptr[rw_idx + 1];
            //std::cout << start << " " << end << std::endl;
            for(int i=start;i<end;i++){
                ierr = fscanf(fp,"%d%f\n",indices +i ,data + i);
            }
        }
    }

    T operator() (int i,int j){
        int start = indptr[i];
        int end = indptr[i + 1];
        T v = 0.;

        for(int k=start;k<end;k++){
            if(indices[k] == j) {
                v = data[k];
                break;
            }
        }

        return v;
    }


}; 
 
