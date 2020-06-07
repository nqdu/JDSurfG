#pragma once
#include"lsmr.hpp"
#include<iostream>

// coordinate format sparse matrix
template<typename T>
class coo_matrix{
    private:
    bool flag; // allocate or not
    int m,n;

    public:
    int *rw,*col;
    int nonzeros;
    T *val;

    // constructor
    coo_matrix(int rows,int cols,int nar)
    {
        m = rows;
        n = cols;
        nonzeros = nar;
        rw = new int [nonzeros]; col = new int [nonzeros];
        val = new T [nonzeros];
        flag = true;
    }

    //contructor
    coo_matrix()
    {
        flag = false;
    }

    // destructor
    ~coo_matrix(){
        if(flag){
            delete [] rw; 
            delete[] col;
            delete[] val;
            flag = false;
        }
    }

    void initialize(int rows,int cols,int nar)
    {
        m = rows;
        n = cols;
        nonzeros = nar;
        rw = new int [nonzeros]; col = new int [nonzeros];
        val = new T [nonzeros];
        flag = true;
    }
    
    // get rows
    int rows(){
        return m;
    }

    // get cols
    int cols(){
        return n;
    }

    void resize(int rows,int cols,int nar)
    {
        m = rows; n = cols;

        // reallocate space
        int *rw_new = new int [nar];
        int *col_new = new int [nar];
        T *value = new T [nar];

        // copy elements
        int nn = nonzeros;
        if(nn> nar) nn = nar;
        for(int i=0;i<nn;i++){
            rw_new[i] = rw[i];
            col_new[i] = col[i];
            value[i] = val[i];
        }
        delete[] val; val = value;
        delete[] rw; rw = rw_new;
        delete[] col; col = col_new;

    }

    void setZeros(){
    // initialize smat
        for(int i=0;i<nonzeros;i++){
            rw[i] = 0;
            col[i] = 0;
            val[i] = 0.0;
        }
    }

    // convert c-style matrix to fortran-style
    void cpp2fortran(bool inverse=false){
        if(!inverse){
            for(int i=0;i<nonzeros;i++){
                rw[i] += 1;
                col[i] +=1;
            }
        }
        else{
            for(int i=0;i<nonzeros;i++){
                rw[i] -= 1;
                col[i] -=1;
            }
        }
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
            for(int k=0;k<nonzeros;k++){
                y[rw[k]] += val[k] * x[col[k]];
            }
        }
        else{
            for(int k=0;k<nonzeros;k++){
                x[col[k]] += val[k] * y[rw[k]];
            }
        }
    }

    // solve A*x =b
    void LsmrSolver(T *b,T *x,LSMRDict<float> &dict){

        // initialize x
        for(int i=0;i<n;i++) x[i] = 0.0;

        // solve by lsmr module
        LSMR(m,n,nonzeros,rw,col,val,b,dict.damp,dict.atol,dict.btol,
            dict.conlim,dict.itnlim,dict.localSize,x,&dict.istop,
            &dict.itn,&dict.anorm,&dict.acond,&dict.rnorm,
            &dict.arnorm,&dict.xnorm);
    }
};