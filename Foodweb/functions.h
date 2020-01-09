#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include "string"
using namespace std;

// GENERAL FUNCTIONS

void create_table(double **table, const int nrow, const int ncol); // create a two dimensions table

void delete_table(double **table, const int nrow); // delete a two dimension table

// PARALLELISATION

int load_parameters_data(const string file_path, const int line_num);

void load_parameters(const string file_path, const int n_simu, double **params, const bool header); // load the file containing the parameters

void setSeed(const int n_simu, int *seeds); // return a vector of integer to create the seeds of the random generator

// OTHER FUNCTIONS

double setAllometric(double z, double M);

double CV(double mean, double SQ, int nPoint); // calculate the coefficient of variation

double SD(double mean, double SQ, int nPoint); // calculate the standard deviation

//INTEGRATION
// System of ODE, function friend of all classes
int func(double t, const double y[], double f[], void *params); // system of ODE

#endif // FUNCTIONS_H_INCLUDED
