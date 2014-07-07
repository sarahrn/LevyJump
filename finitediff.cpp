#include "finitediff.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <random>
using namespace std;

double finitediff::run(){
    //define variables
    double S=exp(10);
    double K=exp(10);
    double t=1;
    int n=10000;
    int m=80;
    double sigma=.2;
    double s=.3;
    double nu=(-.3);

    //choose .353 and 10.229 for X bounds
    double xl=.35305;
    double xu=10.239;
    double deltat=t/n;
    double deltax=(xu-xl)/m;
    double alpha=deltat/pow(deltax,2);

    //choose -.5 and .5 for the truncated jump bounds
    double bl=-.5;
    double bu=.5;

    //translate jump bounds into indices
    int kl=round(bl/deltax-.5);
    int ku=round(bu/deltax+.5);

    //store indices to use for J_(i,j)=sum_-K^K u[i,k+j]*nu_j
    int K1[1][ku-kl+1];
    int K2[1][ku-kl+1];
    for(int j=0;j<(ku-kl+1);j++){
        K1[1][j]=j+1;
        K2[1][j]=j-ku;
    }

    //Let U be computation grid and u be row at U_(i-1)
    //to be used when we increment i in 1...N
    double U[n+1][m+abs(kl)+abs(ku)+1];
    double u[1][m+abs(kl)+abs(ku)+1];

    //Get alphahat and lambdahat
    double alphahat=0;
    double lambdahat=0;
    for(int j=kl;j<=ku;j++){
        double yj=bl+j*deltax;
        double z1=(j+.5)*deltax;
        double z2=(j-.5)*deltax;
        double vj=(-s*erf((nu-z1)/(sqrt(2)*s)))/(2*sqrt(pow(s,2)))-(-s*erf((nu-z2)/(sqrt(2)*s)))/(2*sqrt(pow(s,2)));
        alphahat=alphahat + (exp(yj)-1)*vj;
        lambdahat=lambdahat+vj;
    }

    //Get nu_j
    double vj[ku-kl+1];
    int counter=0;
    for(int j=kl;j<=ku;j++){
        counter=counter+1;
        double z1=(j+.5)*deltax;
        double z2=(j-.5)*deltax;
        vj[counter-1]=(-s*erf((nu-z1)/(sqrt(2)*s)))/(2*sqrt(pow(s,2)))-(-s*erf((nu-z2)/(sqrt(2)*s)))/(2*sqrt(pow(s,2)));
    }

    //set bounds for U
    for(int q=0;q<(n+1);q++){
        U[q][1+abs(kl)]=0;
        U[q][m+abs(kl)]=max((exp(xu)-K),0.0);
    }
    int v=0;
    for(int x=abs(kl);x<(m+abs(kl));x++){
        v=v+1;
        U[1][x+1]=max(((exp(xl+v*deltax)-K)),0.0);
    }

    //initialize u
    for(int j=0;j<(abs(kl)+abs(ku)+m+1);j++){
        u[1][j]=U[1][j];
    }

    //set rows i in 1...N of U
    for(int i=1;i<(n+1);i++){
        for(int j=(1+abs(kl));j<(m+abs(kl));j++){
            //Get J
            double J=0;
            for(int k=0;k<(ku-kl+1);k++){
                J=J+vj[k]*u[1][j+K2[1][k]];
            }
            J=deltat*J;
            //Get D
            double D=u[1][j-1]*((alpha*pow(sigma,2))/2) - u[1][j]*(alpha*pow(sigma,2)) + u[1][j+1]*((alpha*pow(sigma,2))/2)-
                deltat*(pow(sigma,2)/2+alphahat)*(1/deltax)*(u[1][j+1]-u[1][j])-deltat*lambdahat*u[1][j];
            //Set U[i]
            U[i][j]=u[1][j]+J+D;
        }
        //Reset u
        for(int k=0;k<(abs(kl)+abs(ku)+m+1);k++){
            u[1][k]=U[i][k];
        }
    }

    //Get index that corresponds to S=log(X)
    int finalindex=round((log(S)-xl)/deltax)+abs(kl)-1;

    //Print column of U[n] that corresponds to the correct index
    cout << "Price using Finite Difference Scheme: " << U[n][finalindex] << endl;

    return 0;
}

finitediff::finitediff(){
}

