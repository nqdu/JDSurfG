#pragma once
extern "C"{
void delsph(float flat1,float flon1,float flat2,float flon2,float *del);
void empirical_relation(float *vsz,float *vpz,float *rhoz);
void empirical_deriv(float vp,float vs,float *drda,float *dadb);
}

// gaussian noise
double gaussian();

// Gauss-Legendre points and weights
void gauleg(double x1, double x2, double *x, double *w,int n);

// gravity functions
double geocen2geogralat(double hlat);
double geogra2geocenlat(double hlat);
double compute_length(double colatErad, double lonErad, double RE, 
                    double colatSrad, double lonSrad, double RS );
void cal_minDis(double &mindis, double colatrado, double lonrado, double ro, 
            double colatrads1, double lonrads1, double colatrads2, double lonrads2, 
            double rs );