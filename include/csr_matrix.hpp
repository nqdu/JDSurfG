#pragma once
#include<iostream>

template<typename T>
class csr_matrix{
    private:
    int MATRIX_ROW,MATRIX_COL;
    bool ALLOCATE;
    int nonzeros;

    public:
    int * __restrict__ indices,* __restrict__ indptr;
    T * __restrict__ data;
    int nonzeros;

    // constructor
    csr_matrix(int rows,int cols,int nar)
    {
        MATRIX_ROW = rows;
        MATRIX_COL = cols;
        nonzeros = nar;
        indices = new int [nar];
        data = new T [nar];
        indptr = new int [rows + 1];
        ALLOCATE = true;
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
            delete[] index;
            ALLOCATE = false;
        }
    }

    // get rows
    int rows(){
        return MATRIX_ROW;
    }

    int cols(){
        return MATRIX_COL;
    }

    bool is_alloc(){
        return ALLOCATE;
    }

    int get_nonzeros(){
        return nonzeros;
    }

    void initialize(int rows,int cols,int nar)
    {
        MATRIX_ROW = rows;
        MATRIX_COL = cols;
        nonzeros = nar;
        indices = new int [nar];
        data = new T [nar];
        indptr = new int [rows + 1];
        ALLOCATE = true;
    }

    void 

};