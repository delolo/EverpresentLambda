#ifndef RANDOM_H
#define RANDOM_H

#include <vector>
#include <string>

// Current Time Seeder
void InitRandom();

// Integer Seeder
void InitRandom(int seed);

// Return a uniform double in (0,1)
double rnd();

// Get random gaussian with mean mu and s.d. sigma
double rndGaussian(double mu, double sigma);

// Get random standard gaussian
double rndGaussian();

#endif