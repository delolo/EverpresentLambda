/*************************************************************************
 *  Generates an N-step simulation of fluctuating lambda
 *
 *  Created by Michel Buck on 25/07/2014
 *  Copyright (c) 2014 First year. All rights reserved.
 *
 *
 *  Compilation:      g++ main.cpp -o main.o
 *
 *  Execution:        ./main.o [Number of Steps]
 *                    Example : 
 *                    ./main.o 1000
 *
 *  Dependencies:     None
 *
 *  Third party:      mersenne.cpp
 *
 *  User parameters:  N    - number of steps in the simulation
 *
 *************************************************************************/

#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <math.h>
#include <ctime>
#include "../lib/mersenne.cpp"

// THIS IS A RANDOM ORANGE

// SI UNITS AND CONSTANTS
const double PI = M_PI;
double HBAR    = 1.05457173E-34; // m^2 * kg / s
double CLIGHT  = 299792458.0;    // m / s
double GNEWTON = 6.67384E-11;    // m^3 / kg / s^2
double AGEOFUNIVERSE = 4.3E17;   // s
//double LPLANCK = pow(8.0 * PI * GNEWTON * HBAR / pow(CLIGHT, 3.0), 0.5);
double LPLANCK = 1.0;
//double TPLANCK = pow(8.0 * PI * GNEWTON * HBAR / pow(CLIGHT, 5.0), 0.5);
double TPLANCK = 1.0;

// Simulator Class
class Simulator {

private:
    // RNG
    CRandomMersenne *RanGen;

    // number of steps
    int steps;

    // free parameter ell
    double ell;

    // variable vectors
    double* a;          // scale factor
    double* N;          // number of atoms
    double* V;          // volume
    double* y;          // temp variables sum_k^i(dt/a) used to speed up volume calculation
    double* S;          // action
    double* rhomat;     // matter energy density
    double* rhorad;     // radiation energy density
    double* lambda;     // lambda
    double* tau;        // proper time (along isotropic worldlines)
    double* debug;      // array used for debugging;

    // initial conditions and model parameters
    double a0;          // initial scale factor
    double tau0;        // initial time
    double deltatau;    // time increment (if constant)
    double V0;          // initial volume
    double rhomat0;     // initial matter energy density
    double rhorad0;     // initial radiation energy density
    double lambda0;     // initial dark energy density

public:
    // Class constructor
    Simulator(int steps) {
        this->initialize(steps);
    }

    // Class destructor
    ~Simulator() {
        free(a);
        free(N);
        free(V);
        free(S);
        free(rhomat);
        free(rhorad);
        free(lambda);
        free(tau);
        free(debug);
    }

public:
    void runSimulation() {
        for (int i = 0; i < steps; i++) {
            doStep(i);
            printf("%d:\t tau = %E\t a = %E\t rhorad = %E\t rhomat = %E\t lambda = %E\n", i, tau[i], a[i], rhorad[i], rhomat[i], lambda[i]);
        }
    }

    void doStep(int i) {
        // New scale factor
        a[i + 1] = a[i] * (1.0 + sqrt((rhorad[i] + rhomat[i] + lambda[i]) / 3.0) * (tau[i + 1] - tau[i]));

        // New volume (double checked)
        V[i + 1] = 0.0;
        for(int k = 0; k < i + 1; k++) {
            y[k] += (tau[i + 1] - tau[i]) / a[i];
            V[i + 1] += pow(a[k] * y[k], 3.0) * (tau[k + 1] - tau[k]);
        }
        V[i + 1] = 4.0 * PI / 3.0 * V[i + 1];

        // New Cardinality
        N[i + 1] = V[i + 1] / pow(ell, 4.0);

        // New Action
        S[i + 1] = S[i] + rndGaussian() * sqrt(N[i + 1] - N[i]);

        // New lambda
        lambda[i + 1] = S[i + 1] / V[i + 1];

        // New rho
        rhomat[i + 1] = rhomat0 / pow(a[i + 1], 3.0);
        rhorad[i + 1] = rhorad0 / pow(a[i + 1], 4.0);
    }

        void printToFile() {
        FILE *ofp;
        char outputFilename[] = "lambda.txt";
        ofp = fopen(outputFilename, "w");

        if (ofp == NULL) {
          fprintf(stderr, "Can't open output file %s!\n",
                  outputFilename);
          exit(1);
        }

        for(int i = 0; i < steps; i++) {
            fprintf(ofp, "%E\t%E\n", tau[i], lambda[i]);
        }
    }

private:
    void initialize(int steps) {
        // Set number of steps
        this->steps = steps;

        // Set free paremeter ell
        ell = 3.0 * LPLANCK;

        // Set initial values
        //deltatau = (AGEOFUNIVERSE / TPLANCK) / steps;
        this->deltatau = 10.0 * TPLANCK;
        printf("delta-tau = %E\n", deltatau);
        tau0 = TPLANCK;
        a0 = 1.0;
        V0 = 0.0;
        rhomat0 = 0.5;
        rhorad0 = 0.5;
        lambda0 = 0.0;

        // Allocate memory
        a =      (double*) malloc(steps * sizeof(double));
        N =      (double*) malloc(steps * sizeof(double));
        V =      (double*) malloc(steps * sizeof(double));
        y =      (double*) malloc(steps * sizeof(double));
        S =      (double*) malloc(steps * sizeof(double));
        rhomat = (double*) malloc(steps * sizeof(double));
        rhorad = (double*) malloc(steps * sizeof(double));
        lambda = (double*) malloc(steps * sizeof(double));
        tau =    (double*) malloc(steps * sizeof(double));
        debug =  (double*) malloc(steps * sizeof(double));

        // Initialise vectors
        a[0] = a0;
        lambda[0] = lambda0;
        rhomat[0] = rhomat0 / pow(a[0], 3.0);
        rhorad[0] = rhorad0 / pow(a[0], 4.0);
        tau[0] = tau0;
        V[0] = V0;
        N[0] = V0 / pow(ell, 4.0);
        S[0] = 0.0;
        for (int i = 1; i < steps; i++) {
            tau[i] = tau0 + i * deltatau;
            y[i] = 0.0;
        }
        debug[0] = 0.0;

        // Seed RNG
        RanGen = new CRandomMersenne((int) time(0));
    } 

    double rnd() {
        return RanGen->Random();
    }

    double rndGaussian() {
        double phi = 2 * M_PI * rnd();
        double R = sqrt(2 * log(1 / (1 - rnd())));
        return cos(phi) * R;
    }
};

int main(int argc, const char * argv[]) {
    int steps = std::atoi(argv[1]);
    printf("Running simulation for %d steps:\n", steps);
    Simulator* simulator = new Simulator(steps);
    simulator->runSimulation();
    simulator->printToFile();
    return 0;
}
