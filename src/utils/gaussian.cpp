#include<random>
#include<iostream>
using std:: default_random_engine;
using std:: normal_distribution;
using std:: uniform_real_distribution;
using std:: uniform_int_distribution;

//static default_random_engine e(11);
static default_random_engine e(11);
static normal_distribution<double> n(0,1);

/*
    generate random numbers with standard normal distribution
*/
double gaussian()
{
    return n(e);
}
