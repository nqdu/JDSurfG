#include "math.h"

/**
 * \brief Brocher (2005,BSSA) empirical relations, convert vs to vp/rho
 * \param vsz shear wave velocity,km/s
 * \param vpz primary wave velocity,km/s
 * \param rhoz density, g/cc
 */
void 
empirical_relation(float vsz,float &vpz,float &rhoz)
{
    vpz = 0.9409 + 2.0947*vsz - 0.8206*pow(vsz,2)+ 
            0.2683*pow(vsz,3) - 0.0251*pow(vsz,4);
    rhoz = 1.6612 * vpz - 0.4721 * pow(vpz,2) + 
            0.0671 * pow(vpz,3) - 0.0043 * pow(vpz,4) + 
            0.000106 * pow(vpz,5);
}

/**
 * \brief Brocher (2005,BSSA) derivative of empirical relations
 * \param vp/vs primary/secondary wave velocity,km/s
 * \param drda drho / dvp
 * \param dadb dvp / dvs
 */
void 
empirical_deriv(float vp,float vs,float &drda,float &dadb)
{
    drda = 1.6612 - 0.4721*2*vp + 0.0671*3*pow(vp,2) - 
           0.0043*4*pow(vp,3) + 0.000106*5*pow(vp,4);
    dadb = 2.0947 - 0.8206*2*vs + 0.2683*3 * pow(vs,2)
           - 0.0251*4*pow(vs,3);
}
