#include "community.h"
#include "functions.h"

#include <iostream>
#include <fstream>
#include <string>
#include "math.h"
using namespace std;

// INTEGRATION
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_errno.h> // for integration
#include <gsl/gsl_matrix.h> // for integration
#include <gsl/gsl_odeiv.h> // for integration

// CREATION OF THE FOODWEB
void Community::setSystem(double qp, double q, double a, double b){
    // BODY MASSES
    _M[0] = 0; // nutrients
    _M[1] = 0; // detritus
    _M[2] = 1; // primary producer
    _M[3] = 10; // herbivore
    _M[4] = 100; // carnivore

    // METABOLIC RATE
    _x[0] = 0; // nutrients
    _x[1] = 0; // detritus
    _x[2] = setAllometric(qp, _M[2]); // primary producer
    _x[3] = setAllometric(q, _M[3]); // herbivore
    _x[4] = setAllometric(q, _M[4]); // carnivore

    // ATTACK RATE
    _a[0] = 0; // nutrients
    _a[1] = 0; // detritus
    _a[2] = 0; // primary producer
    _a[3] = setAllometric(a, _M[3]); // herbivore
    _a[4] = setAllometric(a, _M[4]); // carnivore

    // HANDLING TIME
    _h[0] = 0; // nutrients
    _h[1] = 0; // detritus
    _h[2] = 0; // primary producer
    _h[3] = pow(b,2) * _M[3] / (6 * 8 * _x[3] * (_M[2] / _M[3]) * _M[2]); // herbivore
    _h[4] = pow(b,2) * _M[4] / (6 * 8 * _x[4] * (_M[2] / _M[3]) * _M[3]); // carnivore

    // ASSIMILATION EFFICIENCY
    _e[0] = 1; // nutrients
    _e[1] = 1; // detritus N
    _e[2] = 1; // primary producer
    _e[3] = 0.45; // herbivore
    _e[4] = 0.85; // carnivore

    // C:N RATIO
    _CN[0] = 1; // nutrients
    _CN[2] = 6.6; // primary producer
    _CN[3] = 5; // herbivore
    _CN[4] = 5; // carnivore
    _CN[1] = _CN[2] * _CN[3] * (1 - _e[3]) / (_CN[3] - _CN[2] * _e[3]); // detritus C:N ratio of herbivore dejections

    // INITIAL POPULATIONS
    _pop[0] = 10; // nutrients
    _pop[1] = 10; // detritus
    _pop[2] = 1; // primary producer
    _pop[3] = 0.5; // herbivore 0.5
    _pop[4] = 0.1; // carnivore 0.1

}
int Community::Dynamic(double t, double tFinal, double tRecord, double tStep, double h, double extinctionThreshold, ofstream &fileChronic){
    //DEFINITIONS FOR GNU LIBRARY
    size_t dimension(5);
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45; // integration method, type of the stepping function
    gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, dimension); // creation an instance of the stepping function
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (h, 0.0); // object keeping the local error on each step within an absolute error of 1e-6 and relative error of 0.0
    gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (dimension); // instance of an evolution function for a system of 1 dimensions
    gsl_odeiv_system sys = {func, NULL, dimension, this}; // the Jacobian is useless with this method : use NULL instead

/////////////////
// INTEGRATION //
/////////////////

    double y[5]; // array with the biomass
    for (int i=0; i<5; i++){
        y[i]=_pop[i]; // initialisation of biomasses
    }
    double nPoint(0); // counter of the number of recorded points
    double tPoint(tRecord); // points when data are recorded

    while (t < tFinal){
        int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, tFinal, &h, y); // integration
        if (status != GSL_SUCCESS)
        break;

        // resources limitation
        for (int i=0; i<2; i++){
            if (y[i]<0){
                y[i]=0; // extinction
            }
        }
        // Extinction
        for (int i=2; i<5; i++){
            if (y[i]<extinctionThreshold && y[i]!=0){
                y[i]=0; // extinction
            }
        }

        ///////////////
        // RECORDING //
        ///////////////

        if (t>tPoint){

            // CHRONICS
            //cout << tPoint << endl;
            fileChronic << t << ";" << _I << ";" << _delta << ";" << _d; // write the values of PP, SP and t
            for (int i=0; i<5; i++){
                fileChronic << ";" << y[i]; // write the values of biomasses and metabolic rates
            }
            for (int i=0; i<6; i++){
                fileChronic << ";" << _flux[i]; // write the values of biomasses and metabolic rates
            }
            fileChronic << ";" << _flux[1]+_flux[2]; // total production
            fileChronic << ";" << _flux[4]+_flux[5]; // total recycling
            fileChronic << ";" << _flux[0]+_flux[1]+_flux[2]; // total production
            fileChronic << ";" << _flux[3]+_flux[4]+_flux[5]; // total recycling
            fileChronic << endl;

            tPoint += tStep ; // next point to be recorded
            nPoint++;
        }
    }

    // END !!!
    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);

    return nPoint;
}

// constructor

Community::Community(
                    double I,
                    double L,
                    double delta,
                    double d,
                    double r,
                    double qp,
                    double q,
                    double B,
                    double K,
                    double FR)
                     :
                    _I(I)
                    ,_L(L)
                    ,_delta(delta)
                    ,_d(d)
                    ,_r(r)
                    ,_qp(qp)
                    ,_q(q)
                    ,_B(B)
                    ,_K(K)
                    ,_FR(FR)
                    {}

// destructor
Community::~Community(){}
