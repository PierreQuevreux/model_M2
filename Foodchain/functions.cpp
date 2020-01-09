#include "functions.h"
#include "community.h"

#include "math.h"
#include <iostream>
#include <fstream>
#include <sstream> // to use istringstream and convert string to double
#include "string"
using namespace std;

// INTEGRATION
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_errno.h> // for integration
#include <gsl/gsl_matrix.h> // for integration
#include <gsl/gsl_odeiv.h> // for integration
#include <gsl/gsl_statistics.h> // for mean and variance

// GENERAL FUNCTIONS
void create_table(double **table, const int nrow, const int ncol){
    for(int i=0; i<nrow; i++){
        table[i] = new double[ncol];
        for(int j=0; j<ncol; j++){
            table[i][j] = 0;
        }
    }
}

void delete_table(double **table, const int nrow){
    for (int i=0; i<nrow; i++){
        delete[] table[i];
        table[i] = NULL;
    }
    delete[] table;
}

int load_parameters(const string file_path, int line_num){
    int out;
    // OPENING OF THE FILE
    ifstream file(file_path.c_str()); // identify the file to open
    if (file){} // control the opening of the file
    else {cout << "ERROR: unable to open file at path" << endl;}
    // DECLARATION OF THE VARIABLEs
    string line; // line of the file
    istringstream iss; // variable for the conversion from string to double

    for(int i=1; i<line_num+1; i++){
        getline(file, line); // read the line of 'file' and store it in 'line'
    }
    iss.str( line.c_str() ); // convert 'value' from string to double
    iss >> out; // write ,the value of the parameter in the table
    return (out);
}

void load_table(const string file_path, const int nrow, double **table, const bool header){
    // OPENING OF THE FILE
    ifstream file(file_path.c_str()); // identify the file to open
    if (file){} // control the opening of the file
    else {cout << "ERROR: unable to open file at path" << endl;}
    // DECLARATION OF THE VARIABLEs
    string line; // line of the file
    string letter; // character in the line
    string value; // value of the parameter (concatenation of the values of 'letter')
    istringstream iss; // variable for the conversion from string to double
    int k(0); // counter of columns
    // HEADER
    if (header){
        getline(file, line); // skip the header
    }
    //

    for(int i=0; i<nrow; i++){
        getline(file, line); // read the line of 'file' and store it in 'line'
        value.clear(); // clear value
        k = 0;
        for (unsigned int j=0; j<line.size(); j++){
            letter = line[j];
            if (letter!=";"){
                value += line[j]; // write the data letter by letter
            }
            else {
                iss.str( value.c_str() ); // convert 'value' from string to double
                iss >> table[i][k];  // write ,the value of the parameter in the table
                iss.clear(); // clear iss
                value.clear(); // clear value
                k++; // next variable
            }
        }
        iss.str( value.c_str() ); // convert 'value' from string to double
        iss >> table[i][k]; // write ,the value of the parameter in the table
        iss.clear(); // clear iss
        value.clear(); // clear value
    }
}

void bifurcation(ofstream &fileBifurcation, const string variables[], const int ntimesteps, const int nstart, const int nparams, const int nvariables, double **TS){
    double table_bif[ntimesteps]; // table with local extrema
    for (int i=0; i<ntimesteps; i++){
        table_bif[i]=-10; // initialisation (values cannot be negative)
    }
    int count_bif(0); // counter for table_bif
    int count_loop(0); // counter for the while loop
    for (int j=nparams+1; j<nparams+1+nvariables; j++){ // +1 because of the time column
        for (int i=nstart+1; i<nstart+ntimesteps-1; i++){
            if( (TS[i][j]<=TS[i-1][j] && TS[i][j]<=TS[i+1][j]) || (TS[i][j]>=TS[i-1][j] && TS[i][j]>=TS[i+1][j]) ){
                while (TS[i][j]!=table_bif[count_loop] && count_loop<count_bif){ // scan the recorded extrema to avoid identical points
                    count_loop++;
                }
                if (count_loop==count_bif){ // this value does not exist yet
                    table_bif[count_bif]=TS[i][j]; // add the value in the table
                    count_bif++;
                }
                count_loop=0; // reset the while loop counter
            }
        }
        for (int i=0; i<count_bif; i++){
            for (int k=1; k<nparams+1; k++){
                fileBifurcation << TS[nstart][k] << ";"; // write parameters values
            }
            fileBifurcation << table_bif[i] << ";" << variables[j] << endl; // write the variable and its name
        }
        for (int i=0; i<ntimesteps; i++){
            table_bif[i]=-10; // reset the table
        }
        count_bif=0; // reset the counter
    }
}

void stat(ofstream &file_mean, ofstream &file_CV, const string variables[], const int ntimesteps, const int nstart, const int nparams, const int nvariables, double **TS){
    double data[ntimesteps]; // array with the time series for one variable
    double mean(0); // mean
    double CV(0); // CV

    // Parameters
    for (int j=1; j<nparams; j++){
        file_mean << TS[nstart][j] << ";"; // write parameters values
        file_CV << TS[nstart][j] << ";"; // write parameters values
    }
    file_mean << TS[nstart][nparams]; // write parameters values
    file_CV << TS[nstart][nparams]; // write parameters values
    // mean and CV
    for (int j=nparams+1; j<nparams+1+nvariables; j++){
        for (int i=nstart; i<nstart+ntimesteps; i++){
            data[i-nstart] = TS[i][j];
        }
        mean = gsl_stats_mean(data, 1, ntimesteps);
        CV = gsl_stats_sd_m(data, 1, ntimesteps, mean);
        if (mean>0){
            CV /= mean;
        }
        else {
            CV = 0;
        }
        file_mean << ";" << mean;
        file_CV << ";" << CV;
    }
    file_mean << endl;
    file_CV << endl;
}

double CV(double mean, double SQ, int nPoint){
    return(pow(SQ/nPoint - pow(mean,2),0.5)/mean);
}

// OTHER FUNCTIONS
double setAllometric(double a, double M){
    return (a*pow(M,-0.25));
}

 // INTEGRATION
 // System of ODE, function friend of all classes
int func(double t, const double y[], double f[], void *params){
    Community *c = (Community *) params;

    // Initialisation of structures
    for(int i=0; i<c->_dim; i++){
        f[i] = 0;
    }
    for(int i=0; i<c->_nbFlux; i++){
        c->_flux[i] = 0;
    }
    for(int i=0; i<c->_diversity; i++){
        c->_GD[i] = 0;
    }
    c->_TR = 0;

// ORGANISMS

    // PREDATION
    double pred[c->_diversity-1]; // vector with functional responses
    pred[0] = c->_a[3] * y[3] * pow(y[2],c->_FR) / (1 + c->_a[3] * c->_h[3] * pow(y[2],c->_FR)); // Herbivore
    pred[1] = c->_a[4] * y[4] * pow(y[3],c->_FR) / (1 + c->_a[4] * c->_h[4] * pow(y[3],c->_FR)); // Carnivore
    pred[2] = c->_a[5] * y[5] * pow(y[4],c->_FR) / (1 + c->_a[5] * c->_h[5] * pow(y[4],c->_FR)); // Top predator

    // RESPIRATION
    double resp[c->_diversity]; // vector with biomass loss due to metabolism and density dependent death
    resp[0] = c->_x[2] * y[2] + c->_B[2] * pow(y[2],2); // Plant
    resp[1] = c->_x[3] * y[3] + c->_B[3] * pow(y[3],2); // Herbivore
    resp[2] = c->_x[4] * y[4] + c->_B[4] * pow(y[4],2); // Carnivore
    resp[3] = c->_x[5] * y[5] + c->_B[5] * pow(y[5],2); // Top predator

    // PRIMARY PRODUCER 2
    double U(0); // nutrient uptake
    U = y[2] * c->_g * y[0]/(c->_K + y[0]);
    f[2] += U // plant growth
            - pred[0] // predation
            - resp[0]; // respiration

    // HERBIVORE 3
    f[3] =  c->_e[3] * pred[0] // predation on plants
            - pred[1] // predation
            - resp[1]; // respiration

    // CARNIVORE 4
    f[4] =  c->_e[4] * pred[1] // predation
            - pred[2] // predation
            - resp[2]; // respiration

    // TOP PREDATOR 5
    f[5] =  c->_e[5] * pred[2] // predation
            - resp[3]; // respiration

    // GROWTH/DEATH RATIO OF SPECIES
    if (y[2]>0){
        c->_GD[0] = U / (pred[0] + resp[0]); // plants
    }
    if (y[3]>0){
        c->_GD[1] = c->_e[3] * pred[0] / (pred[1] + resp[1]); // herbivore
    }
    if (y[4]>0){
        c->_GD[2] = c->_e[4] * pred[1] / (resp[2]); // carnivore
    }
    if (y[5]>0){
        c->_GD[3] = c->_e[5] * pred[2] / (resp[3]); // top predator
    }

// RESOURCES

    // NUTRIENTS 0
    f[0] += c->_d * y[1] // detritus decomposition
            + c->_delta * resp[0]/c->_CN[2] // plant direct recycling
            + c->_delta * resp[1]/c->_CN[3] // herbivore direct recycling
            + c->_delta * resp[2]/c->_CN[4]; // carnivore direct recycling
            + c->_delta * resp[3]/c->_CN[5]; // top predator direct recycling
    c->_TR = f[0]; //  total recycling
    f[0] += c->_I + c->_Irecy - c->_L * y[0] // input + leaching
            - U / c->_CN[2]; // plant uptake

    // DETRITUS 1
    if(c->_d == 0){
        f[1] = 0;
    }
    else{
        f[1] += (1 - c->_delta) * resp[0] / c->_CN[2] // plant wastes
                + (1 - c->_e[3]) * pred[0] / c->_CN[1] + (1 - c->_delta) * resp[1] / c->_CN[3] // non ingested food and wastes of herbivores (beware the stoechiometry !)
                + (1 - c->_e[4]) * pred[1] / c->_CN[3] + (1 - c->_delta) * resp[2] / c->_CN[4] // non ingested food and wastes of carnivores
                + (1 - c->_e[5]) * pred[2] / c->_CN[4] + (1 - c->_delta) * resp[3] / c->_CN[5] // non ingested food and wastes of carnivores
                - c->_L * y[1] - c->_d * y[1]; // detritus N
    }

// PRODUCTION

    // PRIMARY PRODUCTION
    c->_flux[0] = U;

    // SECONDARY PRODUCTION
    c->_flux[1] = c->_e[3] * pred[0];
    c->_flux[2] = c->_e[4] * pred[1];
    c->_flux[3] = c->_e[5] * pred[2];

    // RECYCLING
    //c->_flux[3] = (c->_delta + (1 - c->_delta) * c->_d / (c->_d + c->_L)) * resp[0] / c->_CN[2]; // nutrients recycled by plants
    //c->_flux[4] = (c->_delta + (1 - c->_delta) * c->_d / (c->_d + c->_L)) * resp[1] / c->_CN[3] + (1 - c->_e[3]) * c->_d / (c->_d + c->_L) * pred[0] / c->_CN[1]; // nutrients recycled by herbivores
    //c->_flux[5] = (c->_delta + (1 - c->_delta) * c->_d / (c->_d + c->_L)) * resp[2] / c->_CN[4] + (1 - c->_e[4]) * c->_d / (c->_d + c->_L) * pred[1] / c->_CN[3]; // nutrients recycled by carnivores

    // DIRECT RECYCLING
    c->_flux[4] = c->_delta * resp[0]/c->_CN[2]; // plant direct recycling
    c->_flux[5] = c->_delta * resp[1]/c->_CN[3]; // herbivore direct recycling
    c->_flux[6] = c->_delta * resp[2]/c->_CN[4]; // carnivore direct recycling
    c->_flux[7] = c->_delta * resp[3]/c->_CN[5]; // top predator direct recycling

    // INDIRECT RECYCLING
    c->_flux[8] = c->_d * y[1];

    return GSL_SUCCESS;
}

int funcBis(double t, const double y[], double f[], void *params){
    Community *c = (Community *) params;

    for(int i=0; i<5; i++){
        f[i] = 0;
    }
    for(int i=0; i<7; i++){
        c->_flux[i] = 0;
    }
    for(int i=0; i<3; i++){
        c->_GD[i] = 0;
    }
    c->_TR = 0;

// ORGANISMS

    // PREDATION
    double pred[2]; // vector with
    pred[0] = c->_a[3] * y[3] * pow(y[2],c->_FR) / (1 + c->_a[3] * c->_h[3] * pow(y[2],c->_FR)); // Herbivore
    pred[1] = c->_a[4] * y[4] * pow(y[3],c->_FR) / (1 + c->_a[4] * c->_h[4] * pow(y[3],c->_FR)); // Carnivore

    // RESPIRATION
    double resp[3];
    resp[0] = c->_x[2] * y[2] + c->_B[2] * pow(y[2],2); // Plant
    resp[1] = c->_x[3] * y[3] + c->_B[3] * pow(y[3],2); // Herbivore
    resp[2] = c->_x[4] * y[4] + c->_B[4] * pow(y[4],2); // Carnivore

    // PRIMARY PRODUCER 2
    double U(0); // nutrient uptake
    U = y[2] * c->_g * y[0]/(c->_K + y[0]);
    f[2] += U // plant growth
            - pred[0] // predation
            - resp[0]; // respiration

    // HERBIVORE 3
    f[3] =  c->_e[3] * pred[0] // predation on plants
            - pred[1] // predation
            - resp[1]; // respiration

    // CARNIVORE 4
    f[4] =  c->_e[4] * pred[1] // predation
            - resp[2]; // respiration

    // GROWTH/DEATH RATIO OF SPECIES
    if (y[2]>0){
        c->_GD[0] = U / (pred[0] + resp[0]);
    }
    if (y[3]>0){
        c->_GD[1] = c->_e[3] * pred[0] / (pred[1] + resp[1]);
    }
    if (y[4]>0){
        c->_GD[2] = c->_e[4] * pred[1] / (resp[2]);
    }

// RESOURCES

    // NUTRIENTS 0
    f[0] += c->_d * y[1] // detritus self decomposition
            + c->_delta * resp[0]/c->_CN[2] // plant direct recycling
            + c->_delta * ((1 - c->_e[3]) * pred[0] / c->_CN[1] + resp[1]/c->_CN[3]) // herbivore direct recycling
            + c->_delta * ((1 - c->_e[4]) * pred[1] / c->_CN[3] + resp[2]/c->_CN[4]); // carnivore direct recycling
    c->_TR = f[0]; //  total recycling
    f[0] += c->_I + c->_Irecy - c->_L * y[0] // input + leaching
            - U / c->_CN[2]; // plant uptake

    // DETRITUS 1
    if(c->_d == 0){
        f[1] = 0;
    }
    else{
        f[1] += (1 - c->_delta) * resp[0] / c->_CN[2] // plant wastes
                + (1 - c->_delta) * ((1 - c->_e[3]) * pred[0] / c->_CN[1] + resp[1] / c->_CN[3]) // non ingested food and wastes of herbivores (beware the stoechiometry !)
                + (1 - c->_delta) * ((1 - c->_e[4]) * pred[1] / c->_CN[3] + resp[2] / c->_CN[4]) // non ingested food and wastes of carnivores
                - c->_L * y[1] - c->_d * y[1]; // detritus N
    }

// PRODUCTION

    // PRIMARY PRODUCTION
    c->_flux[0] = U;

    // SECONDARY PRODUCTION
    c->_flux[1] = c->_e[3] * pred[0];
    c->_flux[2] = c->_e[4] * pred[1];

    // RECYCLING
    //c->_flux[3] = (c->_delta + (1 - c->_delta) * c->_d / (c->_d + c->_L)) * resp[0] / c->_CN[2]; // nutrients recycled by plants
    //c->_flux[4] = (c->_delta + (1 - c->_delta) * c->_d / (c->_d + c->_L)) * resp[1] / c->_CN[3] + (1 - c->_e[3]) * c->_d / (c->_d + c->_L) * pred[0] / c->_CN[1]; // nutrients recycled by herbivores
    //c->_flux[5] = (c->_delta + (1 - c->_delta) * c->_d / (c->_d + c->_L)) * resp[2] / c->_CN[4] + (1 - c->_e[4]) * c->_d / (c->_d + c->_L) * pred[1] / c->_CN[3]; // nutrients recycled by carnivores

    // DIRECT RECYCLING
    c->_flux[3] = c->_delta * resp[0]/c->_CN[2]; // plant direct recycling
    c->_flux[4] = c->_delta * ((1 - c->_e[3]) * pred[0] / c->_CN[1] + resp[1]/c->_CN[3]); // herbivore direct recycling
    c->_flux[5] = c->_delta * ((1 - c->_e[4]) * pred[1] / c->_CN[3] + resp[2]/c->_CN[4]); // carnivore direct recycling

    // INDIRECT RECYCLING
    c->_flux[6] = c->_d * y[1];

    return GSL_SUCCESS;
}
