#include "community.h"
#include "species.h"
#include "functions.h"

#include <iostream>
#include <fstream>
#include "string"
#include "math.h"
#include <algorithm> // to use sort
using namespace std;

#include <cstdlib> // to pick random numbers
#include <ctime>  // to pick random numbers

// INTEGRATION
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_errno.h> // for integration
#include <gsl/gsl_matrix.h> // for integration
#include <gsl/gsl_odeiv.h> // for integration
#include <gsl/gsl_randist.h> // for random distributions


// CREATION OF THE FOOD-WEB
void Community::setMass(double Mmin, double Mmax){
    _M = new double[_diversity];

    const gsl_rng_type * T; // type of generator
    gsl_rng * r; // generator
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T); // instance of the generator
    gsl_rng_set(r, time(0)); // change the seed of the generator to have different values after each run

    for (int i=0; i<_diversity; i++){
        _M[i] = pow(10,gsl_rng_uniform(r) * (Mmax - Mmin) + Mmin);
    }
    gsl_rng_free(r);
    sort(_M,_M+_diversity); // sort the elements by ascending order
}
void Community::setObject(){
    _nba=0; // counter to have the number of possible interactions
    for (int i=_nbPrimaryProducer; i<_diversity; i++){
        for (int j=0; j<i; j++){
            if (_M[j]<_M[i]*_b && _M[j]>_M[i]*_b*0.1){ // find all the prays within the feeding niche
                _nba++;
            }
        }
    }
    _pop = new double[_dim+_nba];
        for (int i=0; i<_dim+_nba; i++){
            _pop[i]=0;
        }
    _params = new Species*[_dim];
        for (int i=0; i<_dim; i++){
            _params[i] = NULL;
        }

    _TL = new double[_diversity];
        for (int i=0; i<_nbPrimaryProducer; i++){
            _TL[i]=0;
        }
        for (int i=_nbPrimaryProducer; i<_diversity; i++){
            _TL[i]=0;
        }
    _meanTL = new double[_diversity];
        for (int i=0; i<_nbPrimaryProducer; i++){
            _meanTL[i]=0;
        }
        for (int i=_nbPrimaryProducer; i<_diversity; i++){
            _meanTL[i]=0;
        }
    _tExt = new double[_diversity];
        for (int i=0; i<_diversity; i++){
            _tExt[i]=0;
        }
    _biomass = new double[_diversity];
        for (int i=0; i<_diversity; i++){
            _biomass[i]=0;
        }
    _biomassSQ = new double[_diversity];
        for (int i=0; i<_diversity; i++){
            _biomassSQ[i]=0;
        }
    _CV = new double[_diversity];
        for (int i=0; i<_diversity; i++){
            _CV[i]=0;
        }
}
void Community::setSpecies(double q, double qp, double B, double e[], double CN[], double r, double K, double J, double f, double b, double P0[]){
    for (int i=0; i<_nbPrimaryProducer; i++){
        _params[_nbResource+i] = new Species(_M[i], qp, B, e[1], CN[1], _diversity, _nbResource, _nba, "Plant", r, K, J, f, b); // creation of primary producers
        _pop[_nbResource+i] = P0[0]; // 200
    }
    for (int i=_nbPrimaryProducer; i<_diversity; i++){
        _params[_nbResource+i] = new Species(_M[i], q, B, e[2], CN[2], _diversity, _nbResource, _nba, "Consumer", r, K, J, f, b); // creation of primary producers
        _pop[_nbResource+i] = P0[1]; // 200
    }
}
void Community::setForagingEffort(){
    // DEFINE THE NUMBER OF PREYS AND PREDATORS FOR EACH SPECIES
    for (int i=_nbPrimaryProducer; i<_diversity; i++){
            int nbPrey(0); // prey counter
        for (int j=0; j<_diversity; j++){
            if (_M[j]<_M[i]*_b && _M[j]>_M[i]*_b*0.1){ // find all the preys within the feeding niche
                _params[_nbResource+j]->changeNbPred(_params[_nbResource+j]->getNbPred()+1); // increment the number of predators of the prey
                nbPrey++; // increment the prey counter
            }
        }
        _params[_nbResource+i]->changeNbPrey(nbPrey); // change the number of preys
    }

    // DEFINE THE VECTOR _pred
    for (int i=0; i<_diversity; i++){
        _params[_nbResource+i]->setPred();
    }

    // DEFINE THE INTERACTIONS AND FILL IN THE VECTORS OF POINTERS
    int nTotPrey(0); // counter for the _pop vector
    for (int i=_nbPrimaryProducer; i<_diversity; i++){
        _params[_nbResource+i]->setPrey(); // create the vector of preys' pointers and other values
        int nPrey(0); // counter for the _prey vector
        for (int j=0; j<_diversity; j++){
            if (_M[j]<_M[i]*_b && _M[j]>_M[i]*_b*0.1){ // find all the preys within the feeding niche
                _params[_nbResource+i]->addPrey(j,nPrey,_params[_nbResource+j],nTotPrey); // add the prey j to the diet of i
                _params[_nbResource+i]->setHandlingTime(nPrey,_params[_nbResource+j]); // calculate the handling time of j by i
                _params[_nbResource+j]->addPred(nTotPrey,_nbResource+i,_params[_nbResource+i]); // add i to the predators of j
                _pop[_nbResource+_diversity+nTotPrey]=1/(double)_params[_nbResource+i]->getNbPrey();
                nPrey++; // increment the prey counter
                nTotPrey++; // increment the total prey counter
            }
        }
    }
}
void Community::setResources(double R0[]){
    for(int i=0; i<_nbResource; i++){
        _pop[i] = R0[i];
    }
}
// INTEGRATION
void Community::resetVariables(double I, double delta, double d){
    // INPUT PARAMETERS
    _I = I;
    _delta = delta;
    _d = d;
    // RECORDING STRUCTURES
        _pop[1] = 0; // remove detritus
        for (int i=0; i<_diversity; i++){
            if (_pop[_nbResource+i]>0){
                _tExt[i] = 0;
            }
        }
    // Trophic level
        for (int i=0; i<_diversity; i++){
            _TL[i] = 0;
            _meanTL[i] = 0;
        }
    // Biomasses and coefficient of variation
        for (int i=0; i<_diversity; i++){
            _biomass[i] = 0; // biomass
            _biomassSQ[i] = 0; // biomass square
            _CV[i] = 0; // coefficient of variation
        }
    // Fluxes
        for(int i=0; i<3; i++){
            _flux[i] = 0;
            _sumFlux[i] = 0;
            _sumFluxSQ[i] = 0;
            _fluxCV[i] = 0;
        }
        for(int i=0; i<2; i++){
            _recy[i] = 0;
            _meanRecy[i] = 0;
        }
    // RECORDING VARIABLES
        _connectance = 0;
        _TLmax = 0;
        _meanCV = 0;
        _diversityFinal = 0;
        _nbPrimaryProducerFinal = 0;
}
int Community::Dynamic(double t, double tFinal, double tRecord, double tStep, double h, double extinctionThreshold, double interactionThreshold,bool recyRec){
    //DEFINITIONS FOR GNU LIBRARY
    size_t dimension(_dim+_nba);
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45; // integration method, type of the stepping function
    gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, dimension); // creation an instance of the stepping function
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (h, 0.0); // object keeping the local error on each step within an absolute error of 1e-6 and relative error of 0.0
    gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (dimension); // instance of an evolution function for a system of 1 dimensions
    gsl_odeiv_system sys = {func, NULL, dimension, this}; // the Jacobian is useless with this method : use NULL instead

    double y[_dim+_nba]; // initial value
    for (int i=0; i<_dim+_nba; i++){
        y[i]=_pop[i];
    }

/////////////////
// INTEGRATION //
/////////////////

    _nPoint = 0; // initialisation of the number of recorded points counter
    double tPoint(tRecord); // time of the recorded point
    bool rec(recyRec); // to know if data at t = tRecord have been recorded
    double aTot(0); // for foraging effort scaling
    int index(0); // for foraging effort scaling

    while (t < tFinal){
        int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, tFinal, &h, y); // integration
        if (status != GSL_SUCCESS)
        break;
        // Extinction
        if (y[0]<0){
            y[0]=0;
        }
        if (y[1]<0){
            y[1]=0;
        }
        for (int i=0; i<_diversity; i++){
            if (y[_nbResource+i]<extinctionThreshold && y[_nbResource+i]!=0){
                y[_nbResource+i]=0; // extinction
                _TL[i] = 0; // no more TL
                _meanTL[i] = 0; // no more TL
                _biomass[i] = 0; // no more biomass
                _biomassSQ[i] = 0; // no more sum of square of biomass
                _tExt[i] = t;
                for (int j=0; j<_params[_nbResource+i]->getNbPred(); j++){
                    y[_params[_nbResource+i]->getwIDpred(j)] = 0; // remove this species of all diets
                }
                if (i>=_nbPrimaryProducer){
                    for (int j=0; j<_params[2+i]->getNbPrey(); j++){
                        y[_params[_nbResource+i]->getwIDprey(j)] = 0; // remove the predator interactions
                    }
                }
            }
        }
        // Foraging effort scaling
        for (int i=_nbPrimaryProducer; i<_diversity; i++){
            if (y[_nbResource+i]!=0){
                aTot = 0;
                //index = 0;
                for (int j=0; j<_params[_nbResource+i]->getNbPrey(); j++){
                    index = _params[_nbResource+i]->getwIDprey(j);
                    aTot += y[index];
                }
                for (int j=0; j<_params[_nbResource+i]->getNbPrey(); j++){
                    index = _params[_nbResource+i]->getwIDprey(j);
                    if (aTot!=0){
                        y[index] = y[index]/aTot;
                    }
                }
            }
        }

        ///////////////
        // RECORDING //
        ///////////////

        if (t>=tRecord && rec){
            for (int i=0; i<_dim+_nba; i++){
                _pop[i] = y[i]; // record the densities at t=tRecord
            }
            rec = false;
        }

        if (t>tPoint){
            // Biomasses recoding
            //cout << tPoint << endl;
            for (int i=0; i<_diversity; i++){
                _biomass[i] += y[_nbResource+i]; // add the value of the biomass
                _biomassSQ[i] += pow(y[_nbResource+i],2); // sum of square of biomass
            }

            // TROPHIC LEVEL
            for (int i=0; i<_nbPrimaryProducer; i++){
                if (y[_nbResource+i]>0){
                    _TL[i]=1;
                    _meanTL[i]+=_TL[i];
                }
            }
            for (int i=_nbPrimaryProducer; i<_diversity; i++){
                if (y[_nbResource+i]>0){
                    _TL[i]=1;
                    for (int j=0; j<_params[_nbResource+i]->getNbPrey(); j++){
                        _TL[i] += _TL[_params[_nbResource+i]->getIDprey(j)-_nbResource] * y[_params[_nbResource+i]->getwIDprey(j)];
                    }
                    _meanTL[i]+=_TL[i];
                }
            }

            // PRODUCTION AND RECYCLING
            for(int i=0; i<3; i++){
                _sumFlux[i] += _flux[i];
                _sumFluxSQ[i] += pow(_flux[i],2);
            }

            // RECYCLING
            for(int i=0; i<2; i++){
                _meanRecy[i] += _recy[i];
            }

            tPoint = tPoint + tStep;
            _nPoint++;
        }
    }

    // CONNECTIVITY & CONNECTANCE
    for (int i=_dim; i<_nba; i++){
        if (y[i]!=0){
            _connectance++; // number of existing interactions
        }
    }

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);
    return 0;
}
// OUTPUT
void Community::output(double data[], double species[], double CV[], double TL[]){
    // TROPHIC LEVEL
    for (int i=0; i<_nbPrimaryProducer; i++){
        _meanTL[i] /= _nPoint;
        TL[i] = _meanTL[i];
        if (_meanTL[i]>_TLmax){
            _TLmax = _meanTL[i]; // keep the highest TL
        }
    }
    for (int i=_nbPrimaryProducer; i<_diversity; i++){
        _meanTL[i] /= _nPoint;
        TL[i] = _meanTL[i];
        if (_meanTL[i]>_TLmax){
            _TLmax = _meanTL[i]; // keep the highest TL
        }
    }

    // BIOMASS
    for (int i=0; i<_diversity; i++){
        _biomass[i] /= _nPoint;
        species[i] = _biomass[i]; // RECORDING OF THE AVERAGE BIOMASS
    }

    // REMAINING SPECIES
    for (int i=0; i<_nbPrimaryProducer; i++){
        if (_tExt[i] == 0){
            _diversityFinal++;
            _nbPrimaryProducerFinal++;
        }
    }
    for (int i=_nbPrimaryProducer; i<_diversity; i++){
        if (_tExt[i]==0){
            _diversityFinal++;
        }
    }

    // CONNECTIVITY & CONNECTANCE
    if (_diversityFinal>0){
        _connectance /= pow(_diversityFinal,2);
    }
    else {_connectance = 0;}

    // PRODUCTIVITY
    for(int i=0; i<3; i++){
        _sumFlux[i] /= _nPoint;
        if (_sumFlux[i] > pow(10,-15)){
            _fluxCV[i] = pow(_sumFluxSQ[i]/_nPoint - pow(_sumFlux[i],2),0.5)/_sumFlux[i]; // coefficient of variation of primary production
        }
    }

    // RECYCLING
    for(int i=0; i<2; i++){
        _meanRecy[i] /= _nPoint;
    }

    // COEFFICIENT OF VARIATION OF BIOMASSES
    _meanCV = 0;
    for (int i=0; i<_diversity; i++){
        if (_biomass[i] != 0){
            _CV[i] = pow(_biomassSQ[i]/_nPoint-pow(_biomass[i],2),0.5)/_biomass[i]; // coefficient of variation of species
            CV[i] = _CV[i]; // record in the output vector
            _meanCV += _CV[i];
        }
    }
    if (_diversityFinal != 0){
        _meanCV /= _diversityFinal;
    }


// FILE WITH GENERAL DATA
    data[0] = 0; // number of simulation
    data[1] = 0; // I
    data[2] = 0; // delta
    data[3] = 0; // d
    data[4] = _diversity; // initial diversity
    data[5] = _nbPrimaryProducer; // initial number of primary producers
    data[6] = _diversityFinal; //
    data[7] = _nbPrimaryProducerFinal;
    data[8] = _sumFlux[0]; // primary production
    data[9] = _sumFlux[1]; // secondary production
    data[10] = _sumFlux[2]; // recycling
    data[11] = _fluxCV[0]; // primary production coefficent of variation
    data[12] = _fluxCV[1]; // secondary production coefficent of variation
    data[13] = _fluxCV[2]; // recycling coefficent of variation
    data[14] = _TLmax; // maximal trophic level
    data[15] = _connectance; // connectence
    data[16] = _meanCV; // average coefficient of variation of organisms
    data[17] = _meanRecy[0]; // nutrients recycled by primary producers
    data[18] = _meanRecy[1]; // nutrients recycled by secondary producers
}
// constructor

Community::Community(int diversity
                     ,int nbPrimaryProducer
                     ,int nbResource
                     ,int dim
                     ,double I
                     ,double L
                     ,double b
                     ,double delta
                     ,double d
                     ,double A
                     ,double B
                     ,double FR)
    :_diversity(diversity)
    ,_nbPrimaryProducer(nbPrimaryProducer)
    ,_nbResource(nbResource)
    ,_dim(dim)
    ,_I(I)
    ,_L(L)
    ,_b(b)
    ,_delta(delta)
    ,_d(d)
    ,_A(A)
    ,_B(B)
    ,_FR(FR)
    ,_connectance(0)
    ,_TLmax(0)
    ,_meanCV(0)
    ,_diversityFinal(0)
    ,_nbPrimaryProducerFinal(0){
        for(int i=0; i<3; i++){
            _flux[i] = 0;
            _sumFlux[i] = 0;
            _sumFluxSQ[i] = 0;
            _fluxCV[i] = 0;
        }
        for(int i=0; i<2; i++){
            _recy[i] = 0;
            _meanRecy[i] = 0;
        }
}
// destructor
Community::~Community(){
    delete[] _M;
    delete[] _pop;
    for (int i=_nbResource; i<_dim; i++){
        delete _params[i];
    }
    delete[] _params;
    delete[] _TL;
    delete[] _meanTL;
    delete[] _tExt;
    delete[] _biomass;
    delete[] _biomassSQ;
    delete[] _CV;
}
