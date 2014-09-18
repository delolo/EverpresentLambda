#pragma once

//#include <sys/_types/_time_t.h>
#include <cmath>
#include <ctime>
#include <new>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>

#include "../lib/mersenne.cpp"

using namespace std;

const double PI = M_PI;

CRandomMersenne *RanGen;

void InitRandom() {
    RanGen = new CRandomMersenne((int) time(0));
}

void InitRandom(int seed) {
    RanGen = new CRandomMersenne((int) seed);
}

double rnd() {
    return RanGen->Random();
}
// get random gaussian with mean mu and s.d. sigma
double rndGaussian(double mu, double sigma) {
    double phi = 2 * M_PI * rnd();
    double R = sqrt(2 * log(1 / (1 - rnd()))) * sigma;
    return mu + cos(phi) * R;
}

double sqr(double x) {
    return x * x;
}

double min(double a, double b) {
    if (a < b)
    return a;
    return b;
}
;

template<typename T> string vecToString(const std::vector<T> &vec) {
    ostringstream stream;
    stream << "(" << vec[0];
    for (int i = 1; i < vec.size(); i++) {
        stream << ", " << vec[i];
    }
    stream << ")";
    return stream.str();
}
;


