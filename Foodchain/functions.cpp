#include "functions.h"
#include "community.h"

#include "math.h"
#include <iostream>
#include <fstream>
#include "string"
using namespace std;

// INTEGRATION
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_errno.h> // for integration
#include <gsl/gsl_matrix.h> // for integration
#include <gsl/gsl_odeiv.h> // for integration

double setAllometric(double a, double M){
    return (a*pow(M,-0.25));
}

 // INTEGRATION
 // System of ODE, function friend of all classes
int func(double t, const double y[], double f[], void *params){
    Community *c = (Community *) params;

    for(int i=0; i<5; i++){
        f[i] = 0;
    }
    for(int i=0; i<6; i++){
        c->_flux[i] = 0;
    }
    for(int i=0; i<3; i++){
        c->_GD[i] = 0;
    }

// ORGANISMS

    // PREDATION
    double pred[2]; // vector with
    pred[0] = c->_a[3] * y[3] * pow(y[2],c->_FR) / (1 + c->_a[3] * c->_h[3] * pow(y[2],c->_FR)); // Herbivore
    pred[1] = c->_a[4] * y[4] * pow(y[3],c->_FR) / (1 + c->_a[4] * c->_h[4] * pow(y[3],c->_FR)); // Carnivore

    // RESPIRATION
    double resp[3];
    resp[0] = c->_x[2] * y[2] + c->_B * pow(y[2],2); // Plant
    resp[1] = c->_x[3] * y[3] + c->_B * pow(y[3],2); // Herbivore
    resp[2] = c->_x[4] * y[4] + c->_B * pow(y[4],2); // Carnivore

    // PRIMARY PRODUCER 2
    double U(0); // nutrient uptake
    U = y[2] * c->_r * y[0]/(c->_K + y[0]);
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

// RESOURCES

    // NUTRIENTS 0
    f[0] += c->_I - c->_L * y[0] // input + leaching
            + c->_d * y[1] // detritus self decomposition
            + c->_delta * resp[0]/c->_CN[2] // plant direct recycling
            + c->_delta * resp[1]/c->_CN[3] // herbivore direct recycling
            + c->_delta * resp[2]/c->_CN[4] // carnivore direct recycling
            - U / c->_CN[2]; // plant uptake

    // DETRITUS 1
    f[1] += (1 - c->_delta) * resp[0] / c->_CN[2] // plant wastes
            + ((1 - c->_e[3]) * pred[0] / c->_CN[1] + (1 - c->_delta) * resp[1]) / c->_CN[3] // non ingested food by herbivores (beware the stoechiometry !)
            + ((1 - c->_e[4]) * pred[1] / c->_CN[3] + (1 - c->_delta) * resp[2]) / c->_CN[4] // non ingested food by carnivores
            - c->_L * y[1] - c->_d * y[1]; // detritus N
    if(c->_d == 0){
        f[1] = 0;
    }

// PRODUCTION

    // PRIMARY PRODUCTION
    c->_flux[0] = U;

    // SECONDARY PRODUCTION
    c->_flux[1] = c->_e[3] * pred[0];
    c->_flux[2] = c->_e[4] * pred[1];

    // RECYCLING
    c->_flux[3] = (c->_delta + (1 - c->_delta) * c->_d / (c->_d + c->_L)) * resp[0] / c->_CN[2]; // nutrients recycled by plants
    c->_flux[4] = (c->_delta + (1 - c->_delta) * c->_d / (c->_d + c->_L)) * resp[1] / c->_CN[3] + (1 - c->_e[3]) * c->_d / (c->_d + c->_L) * pred[0] / c->_CN[1]; // nutrients recycled by herbivores
    c->_flux[5] = (c->_delta + (1 - c->_delta) * c->_d / (c->_d + c->_L)) * resp[2] / c->_CN[4] + (1 - c->_e[4]) * c->_d / (c->_d + c->_L) * pred[1] / c->_CN[3]; // nutrients recycled by carnivores

    return GSL_SUCCESS;
}
