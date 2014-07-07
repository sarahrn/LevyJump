#include "fourier.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <random>
using namespace std;

double fourier::run(){
    //set values
    double t=1;
    double x=10;
    double k=10;
    double sigma=.2;
    double maximum=11.0389;
    double n=10000;
    double epsilon=maximum/n;
    double nu=-.3;
    double delta=.3;
    double zeta=-1.4998173;
    double mu = -exp(nu+.5*pow(delta,2))+1.0-.5*pow(sigma,2);

    //calibrate complex numbers
    typedef complex<double> cdouble;

    //we will add to csum during the for loop
    double csum=0;

    //use the riemann sum approximation of the integral
    for(double j=(-maximum);j<=maximum;(j=j+epsilon)){
        cdouble zetaapprox(j,zeta);
        cdouble i(0,1);
        //get h-hat
        cdouble h=-exp(-i*zetaapprox*k+k)/(i*zetaapprox+pow(zetaapprox,2));
        //get new phi with extra term
        //use the closed-form solution to the integral (e^x-1)*nu from -inf to inf
        cdouble phi=mu*i*zetaapprox-.5*pow(zetaapprox,2)*pow(sigma,2)+exp(i*zetaapprox*nu-.5*pow(i*zetaapprox,2)*pow(delta,2))-1.0;
        //increment csum
        csum=csum+1/(2*PI)*real(exp(i*zetaapprox*x+phi*t)*h*epsilon);
    }

    //print results
    cout << "Price of option using Fourier: " << csum << endl;

    return 0;

}

fourier::fourier(){
}


