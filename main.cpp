#include "finitediff.h"
#include "montecarlo.h"
#include "fourier.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <random>
using namespace std;


int main(){
    finitediff f;
    f.run();
    montecarlo m;
    m.run();
    fourier fo;
    fo.run();
    return 0;
}
