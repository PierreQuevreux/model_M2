#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include "string"
using namespace std;

// GENERAL FUNCTIONS
void create_table(double **table, const int nrow, const int ncol); // create a 2 dimensions table

void delete_table(double **table, const int nrow); // delete a 2 dimensions table

int load_parameters(const string file_path, int line_num);

void load_table(const string file_path, const int nrow, double **params, const bool header); // load a text table

void bifurcation(ofstream &fileTimeSeries, const string variable[], const int ntimesteps, const int nstart, const int nparams, const int nvariables, double **TS); // find the extrema for the bifurcation diagram

void stat(ofstream &file_mean, ofstream &file_CV, const string variables[], const int ntimesteps, const int nstart, const int nparams, const int nvariables, double **TS); // compute the mean and the CV

double CV(double mean, double SQ, int nPoint); // returns the coefficient of variation

// OTHER FUNCTIONS
double setAllometric(double a, double M); // set allometric variables

 // INTEGRATION
int func(double t, const double y[], double f[], void *params); // system of ODE
int funcBis(double t, const double y[], double f[], void *params); // system of ODE
#endif // FUNCTIONS_H_INCLUDED
