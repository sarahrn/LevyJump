#include "montecarlo.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <random>
using namespace std;

double montecarlo::run(){

    //calibrate normal generator
    std::mt19937 generator;
    double mean = 0.0;
    double stdev  = 1.0;
    std::normal_distribution<double> normal(mean, stdev);

    //initialize variables
    int N=500000;
    double T=1;
    int m=50;
    double h=T/m;
    double X0=10;
    double sigma=.2;
    double nu=-.3;
    double delta=.3;
    double K=exp(10);
    double mu=-exp(nu+.5*pow(delta,2))+1.0-.5*pow(sigma,2);

    //set total sum (which will be divided at the end to obtain the average) to zero
    //create array to hold discounted call prices
    double sum=0;
    double farray[N];

    //set control variate sum to be averaged at the end and array to hold control variate values
    //for the control variate we use e^{-rt}*Z-Z0=e^{-rt}*Z-E[e^{-rt}Z]=Z-Z0=Z-K (r is 0, K=Z0)
    double cvsum=0;
    double cvarray[N];

    //start loop of simulations
    for(int j=1;j<=N;j++){
        //initialize array X to store paths of X
        double X[m+1];

        //generate path of X=log(Z) using euler scheme
        X[0]=X0;
        for(int k=1;k<(m+1);k++){
          double rn=normal(generator);
          //antithetic sampling: for half of the simulations N, set
          //the random variable to be the negative of itself to
          //introduce negative correlation
          if(j > .5*N){
            rn = -rn;
          }
          X[k]=X[k-1]+mu*h+sigma*sqrt(h)*rn;
        }

        //set Z=exponent of last value of X. This is the final value
        //of the stock Z=e^X
        double Z=exp(X[m]);

        //initialize array tau to hold jump times
        double tauj=0;
        double tau[m+1];
        for(int k=0;k<(m+1);k++){
            tau[k] = 0;
        }

        //generate jump times. if the jump time is greater than T, kill the process
        int counter=0;
        while(tauj < T){
          counter=counter+1;
          double runif = (rand() % 100)/100.0;
          //antithetic sampling
          if(j>.5*N){
            runif=1-runif;
          }
          double R = (-log(runif));
          tauj = tau[counter-1]+R;
          if(tauj < T){
            tau[counter]=tauj;
          }
        }

        //generate lognormal variables and sum them for every jump time
        double P=0;
        counter=counter-1;
        if(counter >= 1){
            double logy[counter];
            for(int k=0;k<counter;k++){
              double yn=normal(generator);
              //antithetic sampling
              if(j>.5*N){
                yn=-yn;
              }
              logy[k]=nu+delta*yn;
            }
            for(int k=0;k<counter;k++){
                P=P+logy[k];
            }
        }

        //add the summed jumps to Z.
        //only need to add to final value of e^X because the jumps are generated
        //independently of the path of X
        Z=Z*exp(P);

        //add the final value of S minus the strike to the total sum
        sum=sum+max((Z-K),0.0);
        farray[j]=max((Z-K),0.0);
        cvsum=cvsum + (Z-K);
        cvarray[j]=Z-K;
    }

    //find b-star
    double sum3=0;
    double sum4=0;
    for(int z=0;z<N;z++){
        sum3=sum3+(cvarray[z]-cvsum/N)*(farray[z]-sum/N);
        sum4=sum4+pow((cvarray[z]-cvsum/N),2);
    }
    double bstar=sum3/sum4;

    //get price = average of max((Z-K),0)-Z-Z0 (Z0=K in the problem)
    //don't need to discount because r is 0 in the problem
    double sum5=0;
    for(int z=0;z<N;z++){
        sum5=sum5+farray[z]-bstar*cvarray[z];
    }
    double price=sum5/N;

    //get monte carlo error = sum of (X-mean(X))^2/(sqrt(N-1)*sqrt(N))
    double sum2=0;
    for(int z=0;z<N;z++){
        sum2=sum2+pow((farray[z]-bstar*cvarray[z]-sum5/N),2);
    }
    double error=sqrt(sum2)/(sqrt(N-1)*sqrt(N)) ;

    //print final values
    cout << "Price of option using Monte Carlo: " << price << endl;
    cout << "Monte Carlo Error: " << error << endl;

    return 0;
}

montecarlo::montecarlo(){
}

