#include <cmath>
#include <limits>

/**
 * @brief compute the derivative of Legendre Polynomial of order n
 * 
 * @param n order 
 * @param x current location
 */
static double dlegendre(int n,double x){
    double dleg{};
    if(n == 0){
        dleg = 0.;
    }
    else if (n==1){
        dleg = 1.;
    }
    else{
        double leg_down1 = x, leg_down2 = 1., leg;
        double dleg_down1 = 1., dleg_down2 = 0.;
        for(int i=2;i<=n;i++){
            leg = (2*i-1)*x*leg_down1/i - (i-1)*leg_down2/i;
            dleg = dleg_down2 + (2*i-1)*leg_down1;
            leg_down2 = leg_down1;
            leg_down1 = leg;
            dleg_down2 = dleg_down1;
            dleg_down1 = dleg;
        }
    }

    return dleg;
}

/**
 * @brief compute the value of Legendre Polynomial of order n
 * 
 * @param n order 
 * @param x current location
 */
static double legendre(int n,double x){
    double leg{};
    if(n == 0){
        leg = 1;
    }
    else if(n==1){
        leg = x;
    } 
    else{
        double leg_down1 = x, leg_down2 = 1.;
        for(int i=2;i<=n;i++){
            leg = (2*i-1)*x*leg_down1/i - (i-1)*leg_down2/i;
            leg_down2 = leg_down1;
            leg_down1 = leg;
        }
    }

    return leg;
}

void gauleg(double x1, double x2, double *x, double *w,int length)
{
    int n = length - 1;
    const double tolerance = 4.0 * std::numeric_limits<double>::epsilon();
    const int nnewton_iter = 100;
    double leg, dleg, delta;
    switch (n)
    {
    case 0:
        x[0] = 0.;
        w[0] = 2.;
        break;
    case 1:
        x[0] = -std::sqrt(1. /3.);
        x[1] = -x[0];
        w[0] = 1; w[1] = 1.;
        break;
    
    default:

        for(int i = 0; i <= (n+1)/2-1;i++){
            x[i] = -std::cos((2.*i+1)/(2.*n+2.) * M_PI);
            for(int j = 1; j <= nnewton_iter; j ++) {
                leg = legendre(n+1,x[i]);
                dleg = dlegendre(n+1,x[i]);
                delta = -leg/dleg;
                x[i] += delta;

                if ( std::abs(delta) <= tolerance * std::abs(x[i]) ) break;
            }

            // utilize symmetric
            x[n-i] = -x[i];
            dleg = dlegendre(n+1,x[i]);
            w[i] = 2. /((1- x[i] * x[i])*dleg * dleg);
            w[n-i] = w[i]; 
        }

        if(n%2 == 0){
            x[n/2] = 0.;
            dleg = dlegendre(n+1,0.);
            w[n/2] = 2. / dleg / dleg;
        }
        break;
    }

    for(int i = 0; i < length; i ++){
        x[i] = 0.5 * (x2 - x1) * x[i] + 0.5 * (x1 + x2);
        w[i] = 0.5 * (x2 - x1) * w[i];
    }
}
