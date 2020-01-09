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
void Community::setSystem(double r, double qp, double q, double J, double B, double a, double b){
    // BODY MASSES
    _M = new double[_dim];
    _M[0] = 0; // nutrients
    _M[1] = 0; // detritus
    //_M[2] = 1; // primary producer
    //_M[3] = 10; // herbivore
    //_M[4] = 100; // carnivore
    //_M[2] = pow(10,-4); // primary producer
    //_M[3] = pow(10,-2.5); // herbivore
    //_M[4] = pow(10,-1); // carnivore
    _M[2] = pow(10,-4); // primary producer -4
    _M[3] = 40*_M[2]; // herbivore
    _M[4] = 40*_M[3]; // carnivore
    _M[5] = 40*_M[4]; // top predator

    // PRIMARY PRODUCER GROWTH RATE
    _g = setAllometric(r, _M[2]);

    // METABOLIC RATE
    _x = new double[_dim];
    _x[0] = 0; // nutrients
    _x[1] = 0; // detritus
    _x[2] = setAllometric(qp, _M[2]); // primary producer
    _x[3] = setAllometric(q, _M[3]); // herbivore
    _x[4] = setAllometric(q, _M[4]); // carnivore
    _x[5] = setAllometric(q, _M[5]); // top predator

    // DENSITY DEPENDENCE MORTALITY
    _B = new double[_dim];
    _B[0] = 0; // nutrients
    _B[1] = 0; // detritus
    _B[2] = setAllometric(B, _M[2]); // primary producer
    _B[3] = setAllometric(B, _M[3]); // herbivore
    _B[4] = setAllometric(B, _M[4]); // carnivore
    _B[5] = setAllometric(B, _M[5]); // top predator

    // ATTACK RATE
    _a = new double[_dim];
    _a[0] = 0; // nutrients
    _a[1] = 0; // detritus
    _a[2] = 0; // primary producer
    _a[3] = setAllometric(a, _M[3]); // herbivore
    _a[4] = setAllometric(a, _M[4]); // carnivore
    _a[5] = setAllometric(a, _M[5]); // top predator

    // HANDLING TIME
    _h = new double[_dim];
    _h[0] = 0; // nutrients
    _h[1] = 0; // detritus
    _h[2] = 0; // primary producer
    _h[3] = pow(b,2) * _M[3] / (6 * J * _x[3] * (b - _M[2] / _M[3]) * _M[2]); // herbivore
    _h[4] = pow(b,2) * _M[4] / (6 * J * _x[4] * (b - _M[3] / _M[4]) * _M[3]); // carnivore
    _h[5] = pow(b,2) * _M[5] / (6 * J * _x[5] * (b - _M[4] / _M[5]) * _M[4]); // top predator

    // ASSIMILATION EFFICIENCY
    _e = new double[_dim];
    _e[0] = 1; // nutrients
    _e[1] = 1; // detritus N
    _e[2] = 1; // primary producer
    _e[3] = 0.45; // herbivore
    _e[4] = 0.85; // carnivore
    _e[5] = 0.85; // top predator

    // C:N RATIO
    _CN = new double[_dim];
    _CN[0] = 1; // nutrients
    _CN[2] = 6.6; // primary producer
    _CN[3] = 5; // herbivore
    _CN[4] = 5; // carnivore
    _CN[5] = 5; // top predator
    _CN[1] = _CN[2] * _CN[3] * (1 - _e[3]) / (_CN[3] - _CN[2] * _e[3]); // detritus C:N ratio of herbivore non ingested food

    // INITIAL POPULATIONS
    _pop = new double[_dim];
    _pop[0] = 10; // nutrients
    _pop[1] = 10; // detritus
    _pop[2] = 1; // primary producer
    _pop[3] = 0.5; // herbivore 0.5
    _pop[4] = 0.1; // carnivore 0.1
    _pop[5] = 0.1; // top predator 0.1

    // OTHER
    _flux = new double[_nbFlux];
    _GD = new double[_diversity];
    _mean = new double[_nvariables];
    _CV = new double[_nvariables];
}

void Community::resetVariables(double delta, double d, bool corr){
    // INPUT PARAMETERS
    _delta = delta;
    _d = d;
    _corr = corr;
    // RECORDING STRUCTURES
    _TR = 0;
    _meanTR = 0;
}

int Community::Dynamic(double t, double tFinal, double tRecord, double tStep, double h, double extinctionThreshold,
                       string parameterValue, ofstream &fileTimeSeries, ofstream &fileMean, ofstream &fileCV){
    //DEFINITIONS FOR GNU LIBRARY
    size_t dimension(_dim);
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45; // integration method, type of the stepping function
    gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, dimension); // creation an instance of the stepping function
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (h, 0.0); // object keeping the local error on each step within an absolute error of 1e-6 and relative error of 0.0
    gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (dimension); // instance of an evolution function for a system of 1 dimensions
    gsl_odeiv_system sys = {func, NULL, dimension, this}; // the Jacobian is useless with this method : use NULL instead

/////////////////
// INTEGRATION //
/////////////////

    double y[_dim]; // array with the biomass
    for (int i=0; i<_dim; i++){
        y[i]=_pop[i]; // initialisation of biomasses
    }
    for (int i=0; i<_nvariables; i++){
        _mean[i]=0;
        _CV[i]=0;
    }
    double nPoint(0); // counter of the number of recorded points
    double tPoint(tRecord); // points at which data are recorded
    int countvar(0); // counter to fill in the _mean and _CV vectors

    while (t < tFinal){
        int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, tFinal, &h, y); // integration
        if (status != GSL_SUCCESS)
        break;

        // resources limitation
        for (int i=0; i<_nbResource; i++){
            if (y[i]<0){
                y[i]=0; // extinction
            }
        }
        // Extinction
        for (int i=_nbResource; i<_dim; i++){
            if (y[i]<extinctionThreshold && y[i]!=0){
                y[i]=0; // extinction
            }
        }

        ///////////////
        // RECORDING //
        ///////////////

        if (t>tPoint){

            // MEAN AND CV
            countvar = 0;
            for (int i=0; i<_dim; i++){
                _mean[i+countvar] += y[i]; // species biomass
                _CV[i+countvar] += pow(y[i],2); // species biomass
            }
            countvar += _dim;
            for (int i=0; i<_diversity; i++){
                _mean[i+countvar] += _GD[i]; // growth:death ratio of species
                _CV[i+countvar] += pow(_GD[i],2); // growth:death ratio of species
            }
            countvar += _diversity;
            for (int i=0; i<_nbFlux; i++){
                _mean[i+countvar] += _flux[i]; // values of production and recycling at species level
                _CV[i+countvar] += pow(_flux[i],2); // values of production and recycling at species level
            }
            countvar += _nbFlux;
            _mean[0+countvar] += _flux[1]+_flux[2]+_flux[3]; // secondary production
            _CV[0+countvar] += pow(_flux[1]+_flux[2]+_flux[3],2); // secondary production
            _mean[1+countvar] += _flux[5]+_flux[6]+_flux[7]; // secondary recycling
            _CV[1+countvar] += pow(_flux[5]+_flux[6]+_flux[7],2); // secondary recycling
            _mean[2+countvar] += _flux[0]+_flux[1]+_flux[2]+_flux[3]; // total production
            _CV[2+countvar] += pow(_flux[0]+_flux[1]+_flux[2]+_flux[3],2); // total production
            _mean[3+countvar] += _flux[4]+_flux[5]+_flux[6]+_flux[7]; // total direct recycling
            _CV[3+countvar] += pow(_flux[4]+_flux[5]+_flux[6]+_flux[7],2); // total direct recycling
            _mean[4+countvar] += _TR; // total recycling
            _CV[4+countvar] += pow(_TR,2); // total recycling

            // CHRONICS
            fileTimeSeries << t << ";" << parameterValue;
            for (int i=0; i<_dim; i++){
                fileTimeSeries << ";" << y[i]; // write the values of biomasses
            }
            for (int i=0; i<_diversity; i++){
                fileTimeSeries << ";" << _GD[i]; // growth:death ratio of species
            }
            for (int i=0; i<_nbFlux; i++){
                fileTimeSeries << ";" << _flux[i]; // write the values of production and recycling at species level
            }
            fileTimeSeries << ";" << _flux[1]+_flux[2]+_flux[3]; // secondary production
            fileTimeSeries << ";" << _flux[5]+_flux[6]+_flux[7]; // secondary recycling
            fileTimeSeries << ";" << _flux[0]+_flux[1]+_flux[2]+_flux[3]; // total production
            fileTimeSeries << ";" << _flux[4]+_flux[5]+_flux[6]+_flux[7]; // total direct recycling
            fileTimeSeries << ";" << _TR; // total recycling
            fileTimeSeries << endl;
            _meanTR += _TR; // average recycling

            tPoint += tStep ; // next point to be recorded
            nPoint++;
        }
    }
    _meanTR /= nPoint;
    _Irecy = _meanTR; // correction for the next simulation

    fileMean << parameterValue;
    fileCV << parameterValue;
    for (int i=0; i<_nvariables; i++){
        _mean[i] /= nPoint; // compute the mean
        _CV[i] = CV(_mean[i], _CV[i], nPoint); // compute the CV
        fileMean << ";" << _mean[i];
        fileCV << ";" << _CV[i];
    }
    fileMean << endl;
    fileCV << endl;

    // END !!!
    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);

    return nPoint;
}

// constructor

Community::Community(
                    int simu_ID,
                    int diversity,
                    int nbResource,
                    int dim,
                    int nbFlux,
                    int nvariables,
                    double I,
                    double L,
                    double delta,
                    double d,
                    double K,
                    double FR)
                     :
                    _simu_ID(simu_ID)
                    ,_diversity(diversity)
                    ,_nbResource(nbResource)
                    ,_dim(dim)
                    ,_nbFlux(nbFlux)
                    ,_nvariables(nvariables)
                    ,_I(I)
                    ,_L(L)
                    ,_delta(delta)
                    ,_d(d)
                    ,_K(K)
                    ,_FR(FR)
                    ,_TR(0)
                    ,_meanTR(0)
                    ,_Irecy(0)
                    ,_corr(false)
                    {}

// destructor
Community::~Community(){
    delete[] _M; // body masses of species
    delete[] _x; // metabolic rate
    delete[] _B; // density dependent mortality rate
    delete[] _a; // attack rate
    delete[] _h; // handling time
    delete[] _e; // assimilation efficiency
    delete[] _CN; // C:N ratio
    delete[] _pop; // initial biomasses
    delete[] _flux; // array with the primary production, herbivore production, carnivore production, primary prod. direct recycling, herbivore direct recycling and carnivore direct recycling and indirect recycling
    delete[] _GD;
    delete[] _mean;
    delete[] _CV;
}
