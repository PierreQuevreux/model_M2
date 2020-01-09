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
#include <gsl/gsl_errno.h> // for integration
#include <gsl/gsl_matrix.h> // for integration
#include <gsl/gsl_odeiv.h> // for integration
#include <gsl/gsl_randist.h> // for random distributions


// CREATION OF THE FOOD-WEB
void Community::setMass(double Mmin, double Mmax, int seed){
    _M = new double[_diversity];

    const gsl_rng_type * T; // type of generator
    gsl_rng * r; // generator
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T); // instance of the generator
    gsl_rng_set(r, seed); // change the seed of the generator to have different values after each run

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
    _meanTL = new double[_diversity];
    _tExt = new double[_diversity];
        for (int i=0; i<_diversity; i++){
            _TL[i]=0;
            _meanTL[i]=0;
            _tExt[i]=0;
        }
    _meanBiomass = new double[_dim];
    _biomassSQ = new double[_dim];
    _biomassCV = new double[_dim];
        for (int i=0; i<_dim; i++){
            _meanBiomass[i]=0;
            _biomassSQ[i]=0;
            _biomassCV[i]=0;
        }
    _recy = new double[_diversity];
    _meanRecy = new double[_diversity];
    _recySQ = new double[_diversity];
    _recyCV = new double[_diversity];
        for (int i=0; i<_diversity; i++){
            _recy[i]=0;
            _meanRecy[i]=0;
            _recySQ[i]=0;
            _recyCV[i]=0;
        }
    _detritus = new double[_diversity];
    _meanDetritus = new double[_diversity];
    _detritusSQ = new double[_diversity];
    _detritusCV = new double[_diversity];
        for (int i=0; i<_diversity; i++){
            _detritus[i]=0;
            _meanDetritus[i]=0;
            _detritusSQ[i]=0;
            _detritusCV[i]=0;
        }
    _flux = new double[_nbFlux];
    _sumFlux = new double[_nbFlux];
    _sumFluxSQ = new double[_nbFlux];
    _fluxCV = new double[_nbFlux];
        for(int i=0; i<_nbFlux; i++){
            _flux[i] = 0;
            _sumFlux[i] = 0;
            _sumFluxSQ[i] = 0;
            _fluxCV[i] = 0;
        }
    _totalBiomass = new double[_nbTotB];
    _meanTotalBiomass = new double[_nbTotB];
    _totalBiomassSQ = new double[_nbTotB];
    _totalBiomassCV = new double[_nbTotB];
        for(int i=0; i<_nbTotB; i++){
            _totalBiomass[i] = 0;
            _meanTotalBiomass[i] = 0;
            _totalBiomassSQ[i] = 0;
            _totalBiomassCV[i] = 0;
        }
    _prod = new double[_nbTotB];
    _meanProd = new double[_nbTotB];
    _prodSQ = new double[_nbTotB];
    _prodCV = new double[_nbTotB];
        for(int i=0; i<_nbTotB; i++){
            _prod[i] = 0;
            _meanProd[i] = 0;
            _prodSQ[i] = 0;
            _prodCV[i] = 0;
        }
    _mort = new double[_nbTotB];
    _meanMort = new double[_nbTotB];
    _mortSQ = new double[_nbTotB];
    _mortCV = new double[_nbTotB];
        for(int i=0; i<_nbTotB; i++){
            _mort[i] = 0;
            _meanMort[i] = 0;
            _mortSQ[i] = 0;
            _mortCV[i] = 0;
        }
}
void Community::setSpecies(double q, double qp, double B, double e[], double CN[], double r, double K, double J, double f, double b, double **pop){
    for (int i=0; i<_nbPrimaryProducer; i++){
        _params[_nbResource+i] = new Species(_M[i], qp, B, e[1], CN[1], _diversity, _nbResource, _nba, "Plant", r, K, J, f, b); // creation of primary producers
        _pop[_nbResource+i] = pop[0][_nbResource+i]; // 10
    }
    for (int i=_nbPrimaryProducer; i<_diversity; i++){
        _params[_nbResource+i] = new Species(_M[i], q, B, e[2], CN[2], _diversity, _nbResource, _nba, "Consumer", r, K, J, f, b); // creation of consumers
        _pop[_nbResource+i] = pop[0][_nbResource+i]; // 5
    }
}
void Community::setForagingEffort(){
    // DEFINE THE NUMBER OF PREY AND PREDATORS FOR EACH SPECIES
    int nbPrey(0); // prey counter
    for (int i=_nbPrimaryProducer; i<_diversity; i++){
        nbPrey = 0; // reset the prey counter
        for (int j=0; j<_diversity; j++){
            if (_M[j]<_M[i]*_b && _M[j]>_M[i]*_b*0.1){ // find all the prey within the feeding niche
                _params[_nbResource+j]->changeNbPred(_params[_nbResource+j]->getNbPred()+1); // increments the number of predators of the prey
                nbPrey++; // increments the prey counter
            }
        }
        _params[_nbResource+i]->changeNbPrey(nbPrey); // set the final number of prey
        _params[_nbResource+i]->setPrey(); // create the vector of prey pointers and other values
    }

    // DEFINE THE VECTOR _pred
    for (int i=0; i<_diversity; i++){
        _params[_nbResource+i]->setPred();
    }

    // DEFINE THE INTERACTIONS AND FILL IN THE VECTORS OF POINTERS
    int nTotPrey(0); // counter for the vector with the foraging efforts (_pop)
    for (int i=_nbPrimaryProducer; i<_diversity; i++){
        nbPrey = 0; // reset the prey counter
        for (int j=0; j<_diversity; j++){
            if (_M[j]<_M[i]*_b && _M[j]>_M[i]*_b*0.1){ // find all the preys within the feeding niche
                _params[_nbResource+i]->addPrey(j,nbPrey,_params[_nbResource+j],nTotPrey); // add the prey j to the diet of i
                _params[_nbResource+i]->setHandlingTime(nbPrey,_params[_nbResource+j]); // calculate the handling time of j by i
                _params[_nbResource+j]->addPred(nTotPrey,_nbResource+i,_params[_nbResource+i]); // add i to the predators of j
                _pop[_nbResource+_diversity+nTotPrey]=1/(double)_params[_nbResource+i]->getNbPrey(); // set the initial foraging effort
                nbPrey++; // increment the prey counter
                nTotPrey++; // increment the total prey counter
            }
        }
    }
}
void Community::setResources(double **pop){
    for(int i=0; i<_nbResource; i++){
        _pop[i] = pop[0][i];
    }
}
// INTEGRATION
void Community::resetVariables(double I, double delta, double d){
    // INPUT PARAMETERS
    _I = I;
    _delta = delta;
    _d = d;
    // RECORDING STRUCTURES
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
        for (int i=0; i<_dim; i++){
            _meanBiomass[i] = 0; // mean biomass
            _biomassSQ[i] = 0; // biomass square
            _biomassCV[i] = 0; // coefficient of variation
        }
    // Fluxes
        for(int i=0; i<_nbFlux; i++){
            _flux[i] = 0;
            _sumFlux[i] = 0;
            _sumFluxSQ[i] = 0;
            _fluxCV[i] = 0;
        }
    // Recycling and coefficient of variation
        for (int i=0; i<_diversity; i++){
            _recy[i] = 0; // recycling
            _meanRecy[i] = 0; // mean recycling
            _recySQ[i] = 0; // recycling square
            _recyCV[i] = 0; // coefficient of variation
            _detritus[i] = 0; // detritus production
            _meanDetritus[i] = 0; // mean detritus production
            _detritusSQ[i] = 0; // detritus production square
            _detritusCV[i] = 0; // coefficient of variation
        }
    // Total biomass, production and mortality
        for(int i=0; i<_nbTotB; i++){
            _totalBiomass[i] = 0;
            _meanTotalBiomass[i] = 0;
            _totalBiomassSQ[i] = 0;
            _totalBiomassCV[i] = 0;
            _prod[i] = 0;
            _meanProd[i] = 0;
            _prodSQ[i] = 0;
            _prodCV[i] = 0;
            _mort[i] = 0;
            _meanMort[i] = 0;
            _mortSQ[i] = 0;
            _mortCV[i] = 0;
        }
    // RECORDING VARIABLES
        _connectance = 0;
        _TLmax = 0;
        _diversityFinal = 0;
        _nbPrimaryProducerFinal = 0;
}
int Community::Dynamic(double t, double tFinal, double tRecord, double tStep, double h, double extinctionThreshold, double interactionThreshold, bool recyRec,
                       string parameterValue, bool time_series_record, ofstream& fileTimeSeriesData, ofstream& fileTimeSeriesBiomass, ofstream& fileTimeSeriesRecy, ofstream& fileTimeSeriesDetritus){
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
    double aTot(0); // sum of foraging efforts of a predator (for scaling)
    int index(0); // position of the foraging effort (for scaling)

    while (t < tFinal){
        int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, tFinal, &h, y); // integration
        if (status != GSL_SUCCESS)
        break;
        // Extinction
        for (int i=0; i<_nbResource; i++){
            if (y[i]<0){
                y[i] = 0; // nutrients or detritus
            }
        }
        for (int i=0; i<_diversity; i++){
            if (y[_nbResource+i]<extinctionThreshold && y[_nbResource+i]!=0){
                y[_nbResource+i] = 0; // extinction
                _TL[i] = 0; // set variables to zero
                _tExt[i] = t;
                _meanTL[i] = 0;
                _meanBiomass[_nbResource+i] = 0;
                _biomassSQ[_nbResource+i] = 0;
                _recy[i] = 0;
                _meanRecy[i] = 0;
                _recySQ[i] = 0;
                _detritus[i] = 0;
                _meanDetritus[i] = 0;
                _detritusSQ[i] = 0;
                for (int j=0; j<_params[_nbResource+i]->getNbPred(); j++){
                    y[_params[_nbResource+i]->getwIDpred(j)] = 0; // remove this species of all diets
                }
                if (i>=_nbPrimaryProducer){
                    for (int j=0; j<_params[_nbResource+i]->getNbPrey(); j++){
                        y[_params[_nbResource+i]->getwIDprey(j)] = 0; // remove the predator interactions
                    }
                }
            }
        }
        // Foraging effort scaling
        for (int i=_nbPrimaryProducer; i<_diversity; i++){
            if (y[_nbResource+i]!=0){
                // calculate the sum of foraging efforts
                aTot = 0;
                for (int j=0; j<_params[_nbResource+i]->getNbPrey(); j++){
                    index = _params[_nbResource+i]->getwIDprey(j);
                    aTot += y[index];
                }
                // rescale all foraging efforts to have values between 0 and 1
                if (aTot!=0){
                    for (int j=0; j<_params[_nbResource+i]->getNbPrey(); j++){
                        index = _params[_nbResource+i]->getwIDprey(j);
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
                _pop[i] = y[i]; // record the densities at t=tRecord for the SC model
            }
            rec = false;
        }

        if (t>tPoint){
            // BIOMASS
            //cout << tPoint << endl;
            for (int i=0; i<_dim; i++){
                _meanBiomass[i] += y[i]; // add the value of the biomass
                _biomassSQ[i] += pow(y[i],2); // sum of square of biomass
            }

            // TOTAL BIOMASS
            for (int i=0; i<_nbTotB; i++){
                _totalBiomass[i] = 0; // reset the total biomass
            }
            for (int i=_nbResource; i<_nbResource+_nbPrimaryProducer; i++){
                _totalBiomass[1] += y[i]; // add the value of the biomass of species i to the total biomass
            }
            for (int i=_nbResource+_nbPrimaryProducer; i<_dim; i++){
                _totalBiomass[2] += y[i]; // add the value of the biomass of species i to the total biomass
            }
            _totalBiomass[0] = _totalBiomass[1] + _totalBiomass[2];
            for (int i=0; i<_nbTotB; i++){
                _meanTotalBiomass[i] += _totalBiomass[i]; // add the value of the biomass of species i to the total biomass
                _totalBiomassSQ[i] += pow(_totalBiomass[i],2); // sum of square of biomass of species i to the total biomass
                _meanProd[i] += _prod[i]; // add the value of the biomass produced by species i to the total biomass
                _prodSQ[i] += pow(_prod[i],2); // sum of square of biomass produced by species i to the total biomass
                _meanMort[i] += _mort[i]; // add the value of the density dependent mortality of species i to the total biomass
                _mortSQ[i] += pow(_mort[i],2); // sum of square of density dependent mortality of species i to the total biomass
            }

            // TROPHIC LEVEL
            for (int i=0; i<_nbPrimaryProducer; i++){
                if (y[_nbResource+i]>0){
                    _TL[i] = 1;
                    _meanTL[i] += _TL[i];
                }
            }
            for (int i=_nbPrimaryProducer; i<_diversity; i++){
                if (y[_nbResource+i]>0){
                    _TL[i] = 1;
                    for (int j=0; j<_params[_nbResource+i]->getNbPrey(); j++){
                        _TL[i] += _TL[_params[_nbResource+i]->getIDprey(j)-_nbResource] * y[_params[_nbResource+i]->getwIDprey(j)];
                    }
                    _meanTL[i] += _TL[i];
                }
            }

            // PRODUCTION AND RECYCLING
            for(int i=0; i<_nbFlux; i++){
                _sumFlux[i] += _flux[i];
                _sumFluxSQ[i] += pow(_flux[i],2);
            }

            // RECYCLING
            for (int i=0; i<_diversity; i++){
                _meanRecy[i] += _recy[i]; // add the value of the recycled nutrients
                _recySQ[i] += pow(_recy[i],2); // sum of square of recycled nutrients
                _meanDetritus[i] += _detritus[i]; // add the value of the produced detritus
                _detritusSQ[i] += pow(_detritus[i],2); // sum of square of produced detritus

            }

//            if (time_series_record){
//                // DATA
//                fileTimeSeriesData << parameterValue << ";" << t;
//                for (int i=0; i<_nbFlux; i++){
//                    fileTimeSeriesData << ";" << _flux[i];
//                }
//                for (int i=0; i<_nbTotB; i++){
//                    fileTimeSeriesData << ";" << _totalBiomass[i];
//                }
//                fileTimeSeriesData << endl;
//                // BIOMASS
//                fileTimeSeriesBiomass << parameterValue << ";" << t;
//                for (int i=0; i<_dim; i++){
//                    fileTimeSeriesBiomass << ";" << y[i];
//                }
//                fileTimeSeriesBiomass << endl;
//                // RECY
//                fileTimeSeriesRecy << parameterValue << ";" << t;
//                for (int i=0; i<_diversity; i++){
//                    fileTimeSeriesRecy << ";" << _recy[i];
//                }
//                fileTimeSeriesRecy << endl;
//                // DETRITUS
//                fileTimeSeriesDetritus << parameterValue << ";" << t;
//                for (int i=0; i<_diversity; i++){
//                    fileTimeSeriesDetritus << ";" << _detritus[i];
//                }
//                fileTimeSeriesDetritus << endl;
//            }

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
void Community::output(double **data, double M[], double **biomass, double **biomassCV, double **recy, double **recyCV, double **detritus, double **detritusCV, double **TL, double **tExt, int run){
    // BODY MASS
    for (int i=0; i<_diversity; i++){
        M[i] = _M[i];
        tExt[i][run] = _tExt[i];
    }

    // BIOMASS
    for (int i=0; i<_dim; i++){
        _meanBiomass[i] /= _nPoint;
        biomass[i][run] = _meanBiomass[i]; // RECORDING OF THE AVERAGE BIOMASS
    }

    // COEFFICIENT OF VARIATION OF BIOMASSES
    for (int i=0; i<_dim; i++){
        if (_meanBiomass[i] != 0){
            _biomassCV[i] = CV(_meanBiomass[i],_biomassSQ[i],_nPoint); // coefficient of variation of species
            biomassCV[i][run] = _biomassCV[i]; // record in the output vector
        }
        else {biomassCV[i][run] = std::numeric_limits<double>::quiet_NaN();} // returns a NaN if the mean is equal to zero
    }

    // RECYCLING AND DETRITUS PRODUCTION
    for (int i=0; i<_diversity; i++){
        _meanRecy[i] /= _nPoint; // average quantity of recycled nutrients
        recy[i][run] = _meanRecy[i];
        _meanDetritus[i] /= _nPoint; // average quantity of produced detritus
        detritus[i][run] = _meanDetritus[i];
        if (_meanRecy[i] != 0){
            _recyCV[i] = CV(_meanRecy[i],_recySQ[i],_nPoint); // coefficient of variation of recycling
            recyCV[i][run] = _recyCV[i]; // record in the output vector
        }
        else {recyCV[i][run] = std::numeric_limits<double>::quiet_NaN();} // returns a NaN if the mean is equal to zero
        if (_meanDetritus[i] != 0){
            _detritusCV[i] = CV(_meanDetritus[i],_detritusSQ[i],_nPoint); // coefficient of variation of produced detritus
            detritusCV[i][run] = _detritusCV[i]; // record in the output vector
        }
        else {detritusCV[i][run] = std::numeric_limits<double>::quiet_NaN();} // returns a NaN if the mean is equal to zero
    }

    // TROPHIC LEVEL
    for (int i=0; i<_diversity; i++){
        _meanTL[i] /= _nPoint;
        TL[i][run] = _meanTL[i];
        if (_meanTL[i]>_TLmax){
            _TLmax = _meanTL[i]; // keep the highest TL
        }
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

    // TOTAL BIOMASS
    for (int i=0; i<_nbTotB; i++){
        _meanTotalBiomass[i] /= _nPoint;
        _meanProd[i] /= _nPoint;
        _meanMort[i] /= _nPoint;
        if (_meanTotalBiomass[i] != 0){
            _totalBiomassCV[i] = CV(_meanTotalBiomass[i],_totalBiomassSQ[i],_nPoint); // coefficient of variation of aggregated species biomass
        }
        else {_totalBiomassCV[i] = std::numeric_limits<double>::quiet_NaN();} // returns a NaN if the mean is equal to zero
        if (_meanProd[i] != 0){
            _prodCV[i] = CV(_meanProd[i],_prodSQ[i],_nPoint); // coefficient of variation of biomass production
        }
        else {_prodCV[i] = std::numeric_limits<double>::quiet_NaN();} // returns a NaN if the mean is equal to zero
        if (_meanMort[i] != 0){
            _mortCV[i] = CV(_meanMort[i],_mortSQ[i],_nPoint); // coefficient of variation of density dependent mortality production
        }
        else {_mortCV[i] = std::numeric_limits<double>::quiet_NaN();} // returns a NaN if the mean is equal to zero
    }

    // CONNECTIVITY & CONNECTANCE
    if (_diversityFinal>0){
        _connectance /= pow(_diversityFinal,2);
    }
    else {_connectance = 0;}

    // GLOBAL RECYCLING
    for(int i=0; i<_nbFlux; i++){
        _sumFlux[i] /= _nPoint;
        if (_sumFlux[i] > pow(10,-15)){
            _fluxCV[i] = CV(_sumFlux[i],_sumFluxSQ[i],_nPoint); // coefficient of variation of flux
        }
        else {_fluxCV[i] = std::numeric_limits<double>::quiet_NaN();} // returns a NaN if the mean is equal to zero
    }

// FILE WITH GENERAL DATA
    data[0][run] = _diversity; // initial diversity (NbSpeciesInit)
    data[1][run] = _nbPrimaryProducer; // initial number of primary producers (NbPPInit)
    data[2][run] = _diversityFinal; // final number of species (NbSpeciesFinal)
    data[3][run] = _nbPrimaryProducerFinal; // final number of primary producers (NbPPFinal)
    data[4][run] = _TLmax; // maximal trophic level (TLmax)
    data[5][run] = _connectance; // connectence
    //
    data[6][run] = _meanTotalBiomass[0]; // average total biomass (totalBiomass)
    data[7][run] = _meanTotalBiomass[1]; // average total biomass of PP (PPbiomass)
    data[8][run] = _meanTotalBiomass[2]; // average total biomass of SP (SPbiomass)
    data[9][run] = _meanProd[0]; // average total biomass production (totalProd)
    data[10][run] = _meanProd[1]; // average total biomass production of PP (PPprod)
    data[11][run] = _meanProd[2]; // average total biomass production of SP (SPprod)
    //
    data[12][run] = _sumFlux[0]; // recycling (Irecy)
    data[13][run] = _sumFlux[1]; // nutrients directly recycled by primary producers (RecyPP)
    data[14][run] = _sumFlux[2]; // nutrients directly recycled by secondary producers (RecySP)
    data[15][run] = _sumFlux[3]; // nutrients directly recycled (RecyDir)
    data[16][run] = _sumFlux[4]; // nutrients indirectly recycled (decomposition of detritus) (RecyInd)
    data[17][run] = _meanMort[0]; // average density mortality (totalMort)
    data[18][run] = _meanMort[1]; // average density mortality of PP (totalMort)
    data[19][run] = _meanMort[2]; // average density mortality of SP (totalMort)
    //
    data[20][run] = _totalBiomassCV[0]; // total biomass coefficient of variation (totalBiomassCV)
    data[21][run] = _totalBiomassCV[1]; // total biomass coefficient of variation of PP (PPbiomassCV)
    data[22][run] = _totalBiomassCV[2]; // total biomass coefficient of variation of SP (SPbiomassCV)
    data[23][run] = _prodCV[0]; // total biomass production coefficient of variation (totalProdCV)
    data[24][run] = _prodCV[1]; // total biomass production coefficient of variation of PP (PPprodCV)
    data[25][run] = _prodCV[2]; // total biomass production coefficient of variation of SP (SPprodCV)
    //
    data[26][run] = _fluxCV[0]; // recycling coefficient of variation (IrecyCV)
    data[27][run] = _fluxCV[1]; // nutrients directly recycled by primary producers coefficient of variation (RecyPPcv)
    data[28][run] = _fluxCV[2]; // nutrients directly recycled by secondary producers coefficient of variation (RecySPcv)
    data[29][run] = _fluxCV[3]; // nutrients directly recycled coefficient of variation (RecyDirCV)
    data[30][run] = _fluxCV[4]; // nutrients indirectly recycled coefficient of variation (decomposition of detritus) (RecyIndCV)
    data[31][run] = _mortCV[0]; // average density mortality coefficient of variation (totalMortCV)
    data[32][run] = _mortCV[1]; // average density mortality of PP coefficient of variation (totalMortCV)
    data[33][run] = _mortCV[2]; // average density mortality of SP coefficient of variation (totalMortCV)
}
// constructor

Community::Community(int diversity
                     ,int nbPrimaryProducer
                     ,int nbResource
                     ,int dim
                     ,int nbFlux
                     ,int nbTotB
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
    ,_nbFlux(nbFlux)
    ,_nbTotB(nbTotB)
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
    ,_diversityFinal(0)
    ,_nbPrimaryProducerFinal(0){}

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
    delete[] _meanBiomass;
    delete[] _biomassSQ;
    delete[] _biomassCV;
    delete[] _recy;
    delete[] _meanRecy;
    delete[] _recySQ;
    delete[] _recyCV;
    delete[] _detritus;
    delete[] _meanDetritus;
    delete[] _detritusSQ;
    delete[] _detritusCV;
    delete[] _flux;
    delete[] _sumFlux;
    delete[] _sumFluxSQ;
    delete[] _fluxCV;
    delete[] _totalBiomass;
    delete[] _meanTotalBiomass;
    delete[] _totalBiomassSQ;
    delete[] _totalBiomassCV;
    delete[] _prod;
    delete[] _meanProd;
    delete[] _prodSQ;
    delete[] _prodCV;
    delete[] _mort;
    delete[] _meanMort;
    delete[] _mortSQ;
    delete[] _mortCV;
}
