#include <cmath>
#include <ctime>
#include <new>
#include <sstream>
#include <iomanip>
#include "../../lib/randomc/mersenne.cpp"

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

double rndGaussian(double mu, double sigma) {
    double phi = 2 * PI * rnd();
    double R = sqrt(2 * log(1 / (1 - rnd()))) * sigma;
    return mu + cos(phi) * R;
}

double rndGaussian() {
        double phi = 2 * PI * rnd();
        double R = sqrt(2 * log(1 / (1 - rnd())));
        return cos(phi) * R;
    }
;
