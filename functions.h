#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

double setAllometric(double z, double M);

//INTEGRATION
// System of ODE, function friend of all classes
int func(double t, const double y[], double f[], void *params); // system of ODE

#endif // FUNCTIONS_H_INCLUDED
