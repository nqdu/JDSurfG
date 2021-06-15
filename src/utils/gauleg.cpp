#include<gsl/gsl_integration.h>

void gauleg(double x1, double x2, double *x, double *w,int n)
{
    gsl_integration_glfixed_table *t;
    t = gsl_integration_glfixed_table_alloc(n);
    
    for(int i=0;i<n;i++){
        gsl_integration_glfixed_point(x1,x2,i,x+i,w+i,t);
    }

    gsl_integration_glfixed_table_free(t);

}