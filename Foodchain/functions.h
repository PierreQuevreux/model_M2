#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include "string"
using namespace std;

double setAllometric(double a, double M);

int func(double t, const double y[], double f[], void *params); // system of ODE

#endif // FUNCTIONS_H_INCLUDED
