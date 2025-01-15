#include <cmath>

const double pi = M_PI;
/*
    Convert geographic latitude to geocentric latitude
    hlat (input) = geographic latitude in radians (north positive)
    glat (output)= geocentric latitude in radians (north positive)
*/
double geogra2geocenlat(double hlat)
{
    double halfpi = 0.5 * pi;
    double x;

    if ((halfpi - fabs(hlat)) >= 0.05) {    //===far away from the pole==
        //x = atan(0.993277*sin(hlat)/cos(hlat));
        x = atan(0.993277*tan(hlat));
    }
    else {
        /* Special formula near pole */
        if (hlat > 0)
        x  = hlat/0.993277 - 0.010632;
        else
        x  = hlat/0.993277 + 0.010632;
    }

    return x;
}

//------convert geocentric latitude back to geographic latitude-------------
//       glatinv (output) = geographic latitude in radians (north positive)
//       hlat (input)= geocentric latitude in radians (north positive)
//---------------------------------------------------------------------
double  geocen2geogralat(double hlat)
{
    double halfpi = 0.5 * pi;
    double glatinv;
    if( (halfpi - fabs(hlat)) >= 0.05){
        //glatinv = atan(sin(hlat)/cos(hlat)/0.993277);
        glatinv = atan(tan(hlat)/0.993277);
    }
    else{   // ------special formula near pole
        if( hlat>0 )
            glatinv = (hlat+0.010632)*0.993277;
        else
            glatinv = (hlat-0.010632)*0.993277;
    }
    return glatinv ;
}

/*
    compute straight line distance between points x_e and x_s points
    the lon,colat,r in rad and km are used as inputs
*/
double compute_length(double colatErad, double lonErad, double RE, 
                    double colatSrad, double lonSrad, double RS )
{
    //double R = 6371.0;
    double PIhalf = 0.5 * pi;
    double Lates, Lones, Res, rep;
    double latErad, latSrad;
    double length;

    //=====Modified 2010/02/15
    latErad = PIhalf - colatErad;
    latSrad = PIhalf - colatSrad;

    Lates   = colatErad - colatSrad;
    Lones   = lonErad - lonSrad;
    Res     = RE - RS;

    rep = Lones * sin( (latErad+latSrad)/2.0 );

    //d2  = Lates*Lates+rep*rep;
    //(*length) = sqrt( Res*Res+d2*RE*RS );
    length = sqrt( Res*Res + (Lates*Lates+rep*rep)*RE*RS );

    return length;
}

/*
    calculate the minmum distance to one grid surface
*/
void cal_minDis(double &mindis, double colatrado, double lonrado, double ro, 
            double colatrads1, double lonrads1, double colatrads2, double lonrads2, 
            double rs )
{
    double length[4];

    /*
        if observation point is outside the same lon/lat
    */
    if( (colatrado<colatrads1||colatrado>colatrads2) && (lonrado<lonrads1||lonrado>lonrads2) ){
        length[0]=compute_length(colatrado, lonrado, ro, colatrads1, lonrads1, rs ); 
        length[1]=compute_length(colatrado, lonrado, ro, colatrads1, lonrads2, rs );
        length[2]=compute_length(colatrado, lonrado, ro, colatrads2, lonrads1, rs );
        length[3]=compute_length(colatrado, lonrado, ro, colatrads2, lonrads2, rs );

        mindis = length[0];
        if( mindis > length[1] ) mindis = length[1];
        if( mindis > length[2] ) mindis = length[2];
        if( mindis > length[3] ) mindis = length[3];
    }
    /*
        if the only longitude of observation point is inside the same lon
    */
    else if((lonrado>=lonrads1&&lonrado<=lonrads2) && (colatrado<colatrads1||colatrado>colatrads2)  ){
        if( fabs(colatrado-colatrads1)<fabs(colatrado-colatrads2) ){
            mindis=compute_length(colatrado, lonrado, ro, colatrads1, lonrado, rs );
        }
        else{
            mindis=compute_length(colatrado, lonrado, ro, colatrads2, lonrado, rs );
        }

    }
    /*
                if the only latitude of observation point is inside the same lat
        */
    else if((colatrado>=colatrads1&&colatrado<=colatrads2) && (lonrado<lonrads1||lonrado>lonrads2) ){
        if( fabs(lonrado-lonrads1)<fabs(lonrado-lonrads2) ){
            mindis=compute_length(colatrado, lonrado, ro, colatrado, lonrads1, rs );
        }
        else{
            mindis = compute_length(colatrado, lonrado, ro, colatrado, lonrads2, rs );
        }
    }
    else{
        mindis = fabs(ro-rs);
    }

}