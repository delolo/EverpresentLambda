#ifndef TOOLS_H
#define TOOLS_H

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

// Square a number
double sqr(double x);

// Find minimum of two integers
double min(double a, double b);

// Convert Vector to String
template<typename T> std::string vecToString(const std::vector<T> &vec);

#endif