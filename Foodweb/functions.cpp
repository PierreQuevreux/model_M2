#include "functions.h"
#include "community.h"
#include "species.h"

#include <iostream>
#include <fstream>
#include <sstream> // to use istringstream and convert string to double

#include "string"
using namespace std;

#include <time.h>
#include "math.h"
#include <gsl/gsl_errno.h> // for integration
#include <gsl/gsl_matrix.h> // for integration
#include <gsl/gsl_odeiv.h> // for integration
#include <gsl/gsl_randist.h> // for random distributions

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

// PARALLELISATION

int load_parameters_data(const string file_path, const int line_num){
    int out;
    // OPENING OF THE FILE
    ifstream file(file_path.c_str()); // identify the file to open
    if (file){} // control the opening of the file
    else {cout << "ERROR: unable to open " << file_path << endl;}
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

void load_parameters(const string file_path, const int n_simu, double **params, const bool header){
    // OPENING OF THE FILE
    ifstream file(file_path.c_str()); // identify the file to open
    if (file){} // control the opening of the file
    else {cout << "ERROR: unable to open " << file_path << endl;}
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

    for(int i=0; i<n_simu; i++){
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
                iss >> params[i][k];  // write ,the value of the parameter in the table
                iss.clear(); // clear iss
                value.clear(); // clear value
                k++; // next variable
            }
        }
        iss.str( value.c_str() ); // convert 'value' from string to double
        iss >> params[i][k]; // write ,the value of the parameter in the table
        iss.clear(); // clear iss
        value.clear(); // clear value
    }
}

void setSeed(const int n_simu, int seeds[]){
    const gsl_rng_type * T; // type of generator
    gsl_rng * r; // generator
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T); // instance of the generator
    gsl_rng_set(r, time(0)); // change the seed of the generator to have different values after each run

    for (int i=0; i<n_simu; i++){
        seeds[i] = gsl_rng_get(r); // draw random integers
    }

    gsl_rng_free(r); // free the memory

}

// OTHER FUNCTIONS

double setAllometric(double z, double M){
    return (z*pow(M,-0.25));
}

double CV(double mean, double SQ, int nPoint){
    return(pow(SQ/nPoint - pow(mean,2),0.5)/mean);
}

double SD(double mean, double SQ, int nPoint){
    return(pow(SQ/nPoint - pow(mean,2),0.5));
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
    for (int i=0; i<c->_nbFlux; i++){
        c->_flux[i] = 0; // Irecy;RecyPP;RecySP;RecyDir;RecyInd
    }
    for (int i=0; i<c->_nbTotB; i++){
        c->_prod[i] = 0; // biomass production
        c->_mort[i] = 0; // density dependent mortality
    }
    for (int i=0; i<c->_diversity; i++){
        c->_recy[i] = 0; // directly recycled nutrients
        c->_detritus[i] = 0; // detritus produced
    }

    // TOTAL OF LOST NUTRIENTS IN BIOMASS DUE TO RESPIRATION AND DENSITY DEPENDENCY
        double loss(0); // total loss of nutrients by all species (expressed as nutrients)
        for (int i=c->_nbResource; i<c->_nbResource+c->_nbPrimaryProducer; i++){
            if(y[i] != 0){
                f[i] -= y[i] * (p[i]->_Q + p[i]->_B*y[i]); // loss for organisms
                c->_mort[1] += p[i]->_B*pow(y[i],2); // density dependent mortality of primary producers
                loss -= f[i] / p[i]->_CN; // lost nutrients (this term must be positive)
                c->_recy[i-c->_nbResource] -= f[i] / p[i]->_CN * c->_delta; // recycling by each primary producer
                c->_detritus[i-c->_nbResource] -= f[i] / p[i]->_CN * (1-c->_delta); // detritus produced by each primary producer
                c->_flux[1] += c->_recy[i-c->_nbResource]; // direct recycling by primary producers
            }
        }
        for (int i=c->_nbResource+c->_nbPrimaryProducer; i<c->_dim; i++){
            if(y[i] != 0){
                f[i] -= y[i] * (p[i]->_Q + p[i]->_B*y[i]); // loss for organisms
                c->_mort[2] += p[i]->_B*pow(y[i],2); // density dependent mortality of consumers
                loss -= f[i] / p[i]->_CN; // lost nutrients (this term must be positive)
                c->_recy[i-c->_nbResource] -= f[i] / p[i]->_CN * c->_delta; // recycling by each secondary producer
                c->_detritus[i-c->_nbResource] -= f[i] / p[i]->_CN * (1-c->_delta); // detritus produced by each primary producer
                c->_flux[2] += c->_recy[i-c->_nbResource]; // recycling by secondary producers
            }
        }
        c->_mort[0] = c->_mort[1] + c->_mort[2]; // total density dependent mortality

    // DETRITUS
        f[1] += (1 - c->_delta) * loss // gain of detritus
        - c->_d * y[1] // decomposition
        - c->_L * y[1]; // leaching

        c->_flux[3] = c->_flux[1] + c->_flux[2]; // direct recycling
        c->_flux[4] = c->_d * y[1]; // indirect recycling

    // NUTRIENTS
        double U(0); // individual plant production
        for (int i=c->_nbResource; i<c->_nbResource+c->_nbPrimaryProducer; i++){
            if(y[i] != 0){
                U = y[i] * p[i]->_r * y[0] / (p[i]->_K + y[0]); // effective growth rate of plants
                f[0] -= U / p[i]->_CN; // nutrient uptake
                f[i] += U; // income for plants
                c->_prod[1] += U; // primary production
            }
        }

        c->_flux[0] += c->_d * y[1] + c->_delta * loss; // recycled nutrients (Irecy)
        f[0] += c->_I // nutrient input
        - c->_L * y[0] // leaching
        + c->_flux[0]; // recycling (Irecy)

    // DENOMINATOR OF FUNCTIONAL RESPONSES OF EACH PREDATOR
        double hPred[c->_dim];
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
        double consumption(0); // eaten biomass of prey i
        double CN(0); // C:N ratio of non ingested biomass (that goes to the detritus)
        for (int k = c->_nbResource + c->_nbPrimaryProducer; k<c->_dim; k++){
            if(y[k] != 0){
                for (int i=0; i<p[k]->_nbPrey; i++){
                    if(y[p[k]->_IDprey[i]] != 0){
                        CN = p[k]->_prey[i]->_CN * p[k]->_CN * (1 - p[k]->_prey[i]->_e) / (p[k]->_CN - p[k]->_prey[i]->_CN * p[k]->_prey[i]->_e); // C:N of the non-ingested biomass
                        consumption = 0;
                        consumption = y[k] * pow(y[p[k]->_IDprey[i]],c->_FR) * y[p[k]->_wIDprey[i]] * p[k]->_a / hPred[k]; // eaten biomass of the prey i
                        f[k] += consumption * p[k]->_prey[i]->_e; // income for the predator
                        c->_prod[2] += consumption * p[k]->_prey[i]->_e; // secondary production
                        f[1] += consumption * (1-p[k]->_prey[i]->_e) / CN; // income for detritus
                        c->_detritus[k-c->_nbResource] += consumption * (1-p[k]->_prey[i]->_e) / CN; // detritus produced by species k
                        f[p[k]->_IDprey[i]] -= consumption; // loss for the prey
                    }
                }

                // ADAPTABILITY
                double growth[p[k]->_nbPrey];
                for (int i=0; i<p[k]->_nbPrey; i++){
                    growth[i]=0;
                    if(y[p[k]->_IDprey[i]] != 0){
                        for (int j=0; j<p[k]->_nbPrey; j++){
                            if(y[p[k]->_IDprey[j]] != 0 && j!=i){
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
        c->_prod[0] = c->_prod[1] + c->_prod[2]; // total biomass production

        if (c->_d==0 && c->_delta==0){
            f[1] = 0; // no change in detritus if no recycling
            c->_flux[0] = 0;
            c->_flux[1] = 0;
            c->_flux[2] = 0;
            c->_flux[3] = 0;
            c->_flux[4] = 0;
        }

    return GSL_SUCCESS;
}
