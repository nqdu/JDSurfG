#pragma once

// spherical distance
void delsph(float colatrad1,float lonrad1,float colatrad2,
             float lonrad2,float &del);

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