#include <omp.h>

void allocate_tasks(int ntasks,int nprocs,int myrank, int &start,int &end)
{
    // allocate jobs to each rank
    int sub_n = ntasks / nprocs;
    int num_larger_procs = ntasks - nprocs * sub_n;
    if (myrank < num_larger_procs){ sub_n = sub_n + 1;
        start = 0 + myrank * sub_n;
    }
    else if (sub_n > 0){ 
        start = 0 + num_larger_procs + myrank * sub_n;
    }
    else { // this process has only zero elements
        start = -1;
        sub_n = 0;
    }
    end = start + sub_n - 1;
}

/**
 * @brief perform y = y + A @ x
 * 
 * @param x shape(col)
 * @param y shape(row)
 * @param indices shape(# of non zeros)
 * @param indptr  shape(row + 1)
 * @param data shape(# of non zeros)
 * @param nprocs # of threads used
 */
static void 
aprod1(const float *x, float * __restrict y,const int *indices,
       const int *indptr, int row,int col,const float *data,int nprocs = 1)
{
    if(nprocs == 1) {
        for(int i = 0;i < row; i++){
        for(int j = indptr[i];j < indptr[i+1];j ++){
            y[i] += data[j] * x[indices[j]];
        }}

        return;
    }

    // backup current omp
    int nprocs_bak = 1;
    #pragma omp parallel 
    {
        nprocs_bak = omp_get_num_threads();
    }
    omp_set_num_threads(nprocs);

    // y += A * x 
    #pragma omp parallel for shared(indptr,indices,x,y,data)
    for(int i = 0;i < row; i++){
        for(int j = indptr[i];j < indptr[i+1];j ++){
            y[i] += data[j] * x[indices[j]];
        }
    }
    omp_set_num_threads(nprocs_bak);
}

/**
 * @brief perform x = x + A.T @ y 
 * 
 * @param x shape(col)
 * @param y shape(row)
 * @param indices shape(# of non zeros)
 * @param indptr  shape(row + 1)
 * @param data shape(# of non zeros)
 * @param nprocs # of threads used
 */
static void 
aprod2(float * __restrict x, const float *y,const int *indices,
       const int *indptr,int row,int col,const float *data,int nprocs = 1)
{
    if(nprocs == 1) {
        for(int i = 0;i < row; i++){
        for(int j = indptr[i];j < indptr[i+1];j ++){
            x[indices[j]] += data[j] * y[i];
        }}

        return;
    }

    // backup current omp
    int nprocs_bak = 1;
    #pragma omp parallel 
    {
        nprocs_bak = omp_get_num_threads();
    }
    omp_set_num_threads(nprocs);

    // y += A * x 
    #pragma omp parallel for shared(indptr,indices,x,y,data)
    for(int i = 0;i < row; i++){
        for(int j = indptr[i];j < indptr[i+1];j ++){
            #pragma omp atomic 
            x[indices[j]] += data[j] * y[i];
        }
    }
    omp_set_num_threads(nprocs_bak);
}

extern "C" {

void myaprod_csr(int mode,float * __restrict x, float * __restrict y,const int *indices,
                const int *indptr,int row,int col, const float *data,int nprocs)
{
    if(mode == 1) {
        aprod1(x,y,indices,indptr,row,col,data,nprocs);
    }
    else {
        aprod2(x,y,indices,indptr,row,col,data,nprocs);
    }
}

}