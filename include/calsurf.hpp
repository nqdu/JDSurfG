extern "C"{

void synthetic(int nx,int ny,int nz,float *vels,float *obst,
            float goxdf,float gozdf,float dvxdf,float dvzdf,
            int kmaxRc,int kmaxRg,int kmaxLc,int kmaxLg,
            double *tRc,double *tRg,double *tLc,double *tLg,int *wavetype,
            int *igrt,int *periods,float *depz,int minthk,float *scxf,
            float *sczf,float *rcxf,float *rczf,int *nrc1,int *nsrcsurf1,
            int kmax,int nsrcsurf,int nrcf);

void CalSurfG(int nx,int ny,int nz,int nonzeros,float *vels,
                  int *rw,int *col,float *val,float *dsurf,
                float goxdf,float gozdf,float dvxdf,float dvzdf,
                int kmaxRc,int kmaxRg,int kmaxLc,int kmaxLg,
                double *tRc,double *tRg,double *tLc,double *tLg,
                int *wavetype,int *igrt,int *periods,float *depz,
                int minthk, float *scxf,float *sczf,float *rcxf,
                float *rczf,int *nrc1,int *nsrcsurf1,int kmax,
                int nsrcsurf,int nrcf,int *nar);
}