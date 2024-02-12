#include <cmath>

/**
 * \brief compute great-circle distance in km between two points
 * \param colatrad1,lonrad1 colat and lon of point1, in rad
 * \param colatrad2,lonrad2 colat and lon of point2, in rad
 * \return del distance in km
 */
extern "C" float delsph(float colatrad1,float lonrad1,float colatrad2,
             float lonrad2)
{
    float R = 6371.0;
    float dlat,dlon,lat1,lat2,a,c;
    float del;
  
    dlat = colatrad2 - colatrad1;
    dlon = lonrad2 - lonrad1;
    lat1 = M_PI * 0.5 - colatrad1;
    lat2 =  M_PI * 0.5 - colatrad2;
    a = sin(dlat * 0.5 ) * sin(dlat * 0.5 ) + sin(dlon * 0.5) 
        * sin(dlon * 0.5) * cos(lat1) * cos(lat2);
    c = 2.0 * atan2(sqrt(a),sqrt(1-a));
    del = R * c;

    return del;
}