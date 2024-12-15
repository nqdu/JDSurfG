#include "shared/csr_matrix.hpp" 
#include <fstream>

extern "C"{
void LSMR(int m, int n, int lenrw, int *rw, int *col, 
        float *val, float *b, float damp, float atol, float btol,
        float conlim, int itnlim, int localSize, float *x, 
        int *istop, int *itn, float *normA, float *condA, 
        float *normr, float *normAr, float *normx);

/**
 * LSMR Solver for Compressed Sparse Row Matrix.       
 * For more details please see the subroutine in src/utils/lsmrModule_csr.f90
 * @param verbose if verbose > 0, print log information on the screen
 * @param num_threads number of threads used in solving linear system
 */
void LSMR_csr(int m, int n, float *val,int *indices,int *indptr, 
            const float *b, float damp, float atol, float btol,
        float conlim, int itnlim, int localSize, float *x, 
        int *istop, int *itn, float *normA, float *condA, 
        float *normr, float *normAr, float *normx,int verbose,
        int num_threads);
}

void csr_matrix:: 
initialize(int rows,int cols,int nar)
{
    // set matrix row/cols 
    MATRIX_ROW = rows;
    MATRIX_COL = cols;
    if(ALLOCATE) {
        delete [] data;
        delete [] indptr;
        delete[] indices;
    }

    nonzeros = nar;
    indices = new int [nar]();
    data = new float [nar]();
    indptr = new int[rows + 1]();
    ALLOCATE = true;

}

/**
 * constructor, set initial values for csr_matrix
*/
csr_matrix::csr_matrix()
{
    this->ALLOCATE = false;
    MATRIX_COL = 0; MATRIX_ROW = 0;
    nonzeros = 0;
    indices = NULL; data = NULL; indptr = NULL;
}

csr_matrix::csr_matrix(int rows,int cols,int nar) 
{
    csr_matrix();
    this -> initialize(rows,cols,nar);
}

csr_matrix::~csr_matrix()
{
    if(ALLOCATE){
        delete [] data;
        delete [] indptr;
        delete[] indices;
        ALLOCATE = false;
        MATRIX_COL = 0; MATRIX_ROW = 0;
        nonzeros = 0;
    }
}

/* 
    matrix-vector multiplication
Parameters:
------------------------------------------
mode : 1.   compute y = y + A * x
    2.   compute x = x + A.T * y
*/
void csr_matrix:: 
aprod(int mode,float *x, float *y) const
{
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

/**
 * @brief convert cpp 0-based to fortran 1-based
 * 
 * @param inverse convert back
 */
void csr_matrix::
cpp2fort(bool inverse)
{
    if(!inverse){
        for(int i=0;i<nonzeros;i++)indices[i] +=1;
        for(int i=0;i<MATRIX_ROW+1;i++) indptr[i] +=1;
    }
    else{
        for(int i=0;i<nonzeros;i++)indices[i] -=1;
        for(int i=0;i<MATRIX_ROW+1;i++) indptr[i] -=1;
    }
}

/**
 * @brief LSMR solver system Ax = b, this function will not change A
 * 
 * @param x unknown variables
 * @param b data
 * @param dict 
 */
void csr_matrix:: 
lsmr_solver(float *x,const float *b,LSMRDict<float> &dict)
{

    // initialize x
    for(int i=0;i<MATRIX_COL;i++) x[i] = 0.0;
    this -> cpp2fort();

    // solve by lsmr module
    int show = 1;
    if(dict.verbose == false) show = 0;

    LSMR_csr(MATRIX_ROW,MATRIX_COL,data,indices,indptr,b,dict.damp,dict.atol,dict.btol,
        dict.conlim,dict.itnlim,dict.localSize,x,&dict.istop,
        &dict.itn,&dict.anorm,&dict.acond,&dict.rnorm,
        &dict.arnorm,&dict.xnorm,show,dict.num_threads);
    this -> cpp2fort(true);
}

/**
 * @brief read binary files to a csr matrix
 * 
 * @param filename name of this bin file
 */
void csr_matrix:: 
read_binary(const std::string &filename)
{
    // scan first to read size of the matrix
    std::ifstream fp(filename,std::ios::binary);
    int nar,nar1,m,n,idx;
    fp.read((char*)&m,sizeof(int));
    fp.read((char*)&n,sizeof(int));

    // read nonzeros/rows in this block
    nar1 = 0;
    while(!fp.read((char*)&idx,sizeof(int))) {
        fp.read((char*)&nar,sizeof(int));
        nar1 += nar;
        
        // read temporay 
        float *tmp = new float[nar * 2];
        fp.read((char*)tmp,sizeof(float)*nar*2);
        delete[] tmp;
    }
    fp.close();

    // allocate space 
    this -> initialize(m,n,nar1);

    // now read data
    fp.open(filename,std::ios::binary);
    fp.read((char*)&m,sizeof(int));
    fp.read((char*)&n,sizeof(int));
    while(!fp.read((char*)&idx,sizeof(int))){
        fp.read((char*)&nar,sizeof(int));
        int start = indptr[idx];
        indptr[idx + 1] = nar + start;
        fp.read((char*)(indices + start),sizeof(int)*nar);
        fp.read((char*)(data + start),sizeof(float)*nar);
    }

    // close file
    fp.close();
}

void csr_matrix:: 
write_binary(const std::string &filename) const
{
    // open file
    std::ofstream fp(filename,std::ios::binary);
    if(!fp.is_open()){
        printf("cannot open file %s\n",filename);
        exit(1);
    } 

    fp.write((char*)&MATRIX_ROW,sizeof(int));
    fp.write((char*)&MATRIX_COL,sizeof(int));
    for(int i = 0; i < MATRIX_ROW; i++){
        int start = indptr[i], end = indptr[i+1];
        int nonzeros = end - start;
        fp.write((char*)&i,sizeof(int));
        fp.write((char*)&nonzeros,sizeof(int));
        fp.write((char*)(indices + start),sizeof(int)*nonzeros);
        fp.write((char*)(data + start),sizeof(int)*nonzeros);
    }

    fp.close();
}

/**
 * merge csr files saved by several procs
 * @param nprocs # of threads used 
 * @param outfile # base filename, such as csr.bin
*/
void merge_csr_files(int nprocs,const std::string &outfile)
{
    std::ifstream fp;
    std::ofstream fpout;
    fpout.open(outfile,std::ios::binary);

    for(int irank = 0; irank < nprocs; irank ++){
        // open filename
        std::string filename = std::string(outfile) + "." + std::to_string(irank);
        fp.open(filename.c_str(),std::ios::binary);

        //get size 
        fp.seekg(0,std::ios::end);
        long size = (size_t)fp.tellg() - sizeof(int) * 2;
        fp.close();

        // read rows/cols 
        fp.open(filename.c_str(),std::ios::binary);
        int m,n;
        fp.read((char*)&m,sizeof(int));
        fp.read((char*)&n,sizeof(int));
        if(irank == 0) {
            fpout.write((char*)&m,sizeof(int));
            fpout.write((char*)&n,sizeof(int));
        }

        // write to fpout
        char *chunk = new char[size + 1];
        fp.read(chunk,size);
        fpout.write(chunk,size);
        fp.close();
        
        // remove this file
        std::remove(filename.c_str());
        delete [] chunk;
    }
    fpout.close();
}

int csr_matrix:: rows() const 
{
    return MATRIX_ROW;
}

int csr_matrix:: cols() const 
{
    return MATRIX_COL;
}

/**
 * Add 2-th order Tikhonov Regularization Matrix to a csr_matrix
 * @param smat     : csr_matrix<float>
 * @param weight   : smooth factor
 * 
 * @example 
 * ----------------------------------------------------
 * For 1-D problem with dimension n = 5, regularization matrix is:
 *             [ 2   0  0  0  0 ]
 *             [ -1  2 -1  0  0 ]
 *             [ 0   -1 2  -1 0 ]
 *             [ 0   0  -1  2 -1]
 *             [ 0   0   0  0  2]
 */
void add_regularization(csr_matrix &smat,float weight,int nx,int ny,int nz)
{
    // get parameters required
    int n = nx * ny * nz; // model dimension
    //int n = (nx-2) * (ny -2) * (nz  - 1); 
    int nar = smat.nonzeros - n * 7; // nonzeros excluding smooth term
    int m = smat.rows() - n; // data dimension

    int count = 0;
    for(int k=0;k<nz;k++){
    for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++){
        if( i==0 || i==nx-1 || j==0 || j==ny-1 || k==0 || k==nz-1){
            
            // and more restrictions to boundary points
            if(nar + 1 > smat.nonzeros){
                printf("please increase sparse ratio!\n");
                exit(1);
            }
            int rwc = count + m;
            smat.data[nar] = 2.0 * weight;
            smat.indices[nar] = k * ny * nx + j * nx + i;
            smat.indptr[rwc + 1] = smat.indptr[rwc] + 1;
            nar += 1;
            count += 1;
            
           continue;
        }
        else{
            if(nar  + 7 > smat.nonzeros){
                printf("please increase sparse ratio!\n");
                exit(1);
            }
            int rwc = count + m;  // current row
            smat.indptr[rwc +1] = smat.indptr[rwc] + 7;
            int clc = k * nx * ny + j * nx + i;// current column
            smat.data[nar] = 6.0 * weight;
            smat.indices[nar] = clc;

            // x direction
            smat.data[nar + 1] = -weight;
            smat.indices[nar + 1] = clc - 1;
            smat.data[nar + 2] = -weight;
            smat.indices[nar + 2] = clc + 1;

            // y direction
            smat.data[nar + 3] = -weight;
            smat.indices[nar + 3] = clc - nx;
            smat.data[nar + 4] = -weight;
            smat.indices[nar + 4] = clc + nx;

            // z direction
            smat.data[nar + 5] = -weight;
            smat.indices[nar + 5] = clc - nx * ny;
            smat.data[nar + 6] = -weight;
            smat.indices[nar + 6] = clc + nx * ny;

            nar += 7;
            count += 1;
        }
    }}}
    smat.nonzeros = nar;
}