#ifndef JDSURFG_GRAVITY_TRANSFORM_H_
#define JDSURFG_GRAVITY_TRANSFORM_H_


// gravity functions
double geocen2geogralat(double hlat);
double geogra2geocenlat(double hlat);
double compute_length(double colatErad, double lonErad, double RE, 
                    double colatSrad, double lonSrad, double RS );
void cal_minDis(double &mindis, double colatrado, double lonrado, double ro, 
            double colatrads1, double lonrads1, double colatrads2, double lonrads2, 
            double rs );

#endif // end JDSURFG_GRAVITY_TRANSFORM_H_
