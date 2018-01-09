#include "functions.h"
#include "community.h"
#include "species.h"

#include "string"
using namespace std;

#include "math.h"
#include <gsl/gsl_errno.h> // for integration
#include <gsl/gsl_matrix.h> // for integration
#include <gsl/gsl_odeiv.h> // for integration


double setAllometric(double z, double M){
    return (z*pow(M,-0.25));
}

//INTEGRATION
// System of ODE, function friend of all classes
int func(double t, const double y[], double f[], void *params){
    Community *c = (Community *) params;
    Species **p = c->_params;

    // Initialization to avoid strange values
    for (int i=0; i<c->_dim+c->_nba; i++){
        f[i]=0;
    }
    for (int i=0; i<c->_diversity; i++){
        c->_popRecy[i] = 0;
    }
    for (int i=0; i<3; i++){
        c->_flux[i] = 0;
    }
    for (int i=0; i<2; i++){
        c->_recy[i] = 0;
    }

    // TOTAL OF LOST NUTRIENTS IN BIOMASS DUE TO RESPIRATION AND DENSITY DEPENDENCY
        double loss(0);
        for (int i=c->_nbResource; i<c->_nbResource+c->_nbPrimaryProducer; i++){
            if(y[i] != 0){
                f[i] -= y[i] * (p[i]->_Q + p[i]->_B*y[i]); // loss for organisms
                loss -= f[i] / p[i]->_CN; // lost nutrients (this term must be positive)
                c->_recy[0] -= f[i] / p[i]->_CN * (c->_delta + (1 - c->_delta) * c->_d / (c->_d + c->_L)); // recycling by primary producers
                c->_popRecy[i - c->_nbResource] -= f[i] / p[i]->_CN * (c->_delta + (1 - c->_delta) * c->_d / (c->_d + c->_L)); // recycling by each species
            }
        }
        for (int i=c->_nbResource+c->_nbPrimaryProducer; i<c->_dim; i++){
            if(y[i] != 0){
                f[i] -= y[i] * (p[i]->_Q + p[i]->_B*y[i]); // loss for organisms
                loss -= f[i] / p[i]->_CN; // lost nutrients (this term must be positive)
                c->_recy[1] -= f[i] / p[i]->_CN * (c->_delta + (1 - c->_delta) * c->_d / (c->_d + c->_L)); // recycling by secondary producers
                c->_popRecy[i - c->_nbResource] -= f[i] / p[i]->_CN * (c->_delta + (1 - c->_delta) * c->_d / (c->_d + c->_L)); // recycling by each species
            }
        }

    // DETRITUS
        f[1] += (1 - c->_delta)*loss // gain of detritus
        - c->_d*y[1] // decomposition
        - c->_L*y[1]; // leaching

    // NUTRIENTS
        double U(0); // individual plant production
        for (int i=c->_nbResource; i<c->_nbResource+c->_nbPrimaryProducer; i++){
            if(y[i] != 0){
                U = y[i]*p[i]->_r*y[0]/(p[i]->_K+y[0]); // effective growth rate of plants
                f[0] -= U / p[i]->_CN; // nutrient uptake
                f[i] += U; // income for plants
                c->_flux[0] += U; // primary production
            }
        }

        c->_flux[2] += c->_d*y[1] + c->_delta*loss; // recycled nutrients
        f[0] += c->_I // nutrient input
        - c->_L*y[0] // leaching
        + c->_flux[2]; // recycling

    // PREDATOR HUNT TIME
    // term (1 + sum(a*f*h*X)), consumption time of all preys by our consumer
        double hPred[c->_dim]; // denominator
        for (int i=c->_nbResource+c->_nbPrimaryProducer; i<c->_dim; i++){
            if(y[i] != 0){
                hPred[i]=0;
                for (int j=0; j<p[i]->_nbPrey; j++){
                    if(y[p[i]->_IDprey[j]] != 0){
                        hPred[i] += y[p[i]->_wIDprey[j]] * p[i]->_a * p[i]->_h[j] * pow(y[p[i]->_IDprey[j]],c->_FR);
                    }
                }
                hPred[i] += 1;
            }
            else{hPred[i]=1;}
        }

    // CONSUMERS
        double consumption(0); // eaten biomass
        double CN(0); // C:N ratio of non ingested biomass (that goes to the detritus)
        for (int k = c->_nbResource + c->_nbPrimaryProducer; k<c->_dim; k++){
            if(y[k] != 0){
                for (int i=0; i<p[k]->_nbPrey; i++){
                    if(y[p[k]->_IDprey[i]] != 0){
                        CN = p[k]->_prey[i]->_CN * p[k]->_CN * (1 - p[k]->_prey[i]->_e) / (p[k]->_CN - p[k]->_prey[i]->_CN * p[k]->_prey[i]->_e); // C:N of the non-ingested biomass
                        consumption = 0;
                        consumption = y[k] * pow(y[p[k]->_IDprey[i]],c->_FR) * y[p[k]->_wIDprey[i]] * p[k]->_a / hPred[k]; // eaten biomass of the prey i
                        f[k] += consumption * p[k]->_prey[i]->_e; // income for the predator
                        c->_flux[1] += consumption * p[k]->_prey[i]->_e; // secondary production
                        f[1] += consumption * (1-p[k]->_prey[i]->_e) / CN; // income for the detritus
                        c->_recy[1] += consumption * (1-p[k]->_prey[i]->_e) / CN * c->_d / (c->_d + c->_L); // recycling by secondary producers
                        c->_popRecy[k - c->_nbResource] += consumption * (1-p[k]->_prey[i]->_e) / CN * c->_d / (c->_d + c->_L); // recycling by each consumer
                        f[p[k]->_IDprey[i]] -= consumption; // loss for the prey
                    }
                }

                // ADAPTABILITY
                double growth[p[k]->_nbPrey];
                for (int i=0; i<p[k]->_nbPrey; i++){
                    growth[i]=0;
                    if(y[p[k]->_IDprey[i]] != 0){
                        for (int j=0; j<p[k]->_nbPrey; j++){
                            if(y[p[k]->_IDprey[j]] != 0){
                                growth[i] += y[p[k]->_wIDprey[j]] * pow(y[p[k]->_IDprey[j]],c->_FR) * p[k]->_a *
                                            (p[k]->_h[j] * p[k]->_prey[i]->_e - p[k]->_h[i] * p[k]->_prey[j]->_e);
                            }
                        }
                        growth[i] += p[k]->_prey[i]->_e;
                        growth[i] *= pow(y[p[k]->_IDprey[i]],c->_FR) * p[k]->_a ;
                    }
                }
                double growthSum(0);
                for (int i=0; i<p[k]->_nbPrey; i++){
                    growthSum += y[p[k]->_wIDprey[i]] * growth[i];
                }
                for (int i=0; i<p[k]->_nbPrey; i++){
                    if(y[p[k]->_wIDprey[i]] != 0){
                        f[p[k]->_wIDprey[i]] = y[p[k]->_wIDprey[i]] * c->_A / (pow(hPred[k],2))
                                                * (growth[i] - growthSum);
                    }
                }
            }
        }

        if (c->_d==0 && c->_delta==0){
            f[1] = 0; // no change in detritus if no recycling
            c->_recy[0] = 0;
            c->_recy[1] = 0;
        }

    return GSL_SUCCESS;
}
