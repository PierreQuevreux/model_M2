#include "community.h"
#include "species.h"
#include "functions.h"

#include <iostream>
#include "math.h"

using namespace std;

int main()
{
string path("data/");

// VARIABLES

// PARAMETERS
// Ecosystem
int diversity(50); // initial number of species
int nbPrimaryProducer(5); // initial number of primary producers
int nbResource(2); // number of resources
int dim(nbResource+diversity); // dimension of the system
double I(40); // nutrient input
double L(0.2); // nutrient leaching rate (arbitrary)
double delta(0.2); // fraction of mineralized excreted nutrients
double d(0.2); // intrinsic mineralization rate of detritus
// Initialisation
double R0[2]{10,10}; // initial density of nutrients and detritus
double P0[2]{10,5}; // initial density of primary producers and predators
// Predator/prey interactions
double b(0.05); //prey/predator body mass ratio limit
double A(0.01); // adaptability
// allometric scaling constants
double r(0.87);  // scaling constant of maximal growth rate
double q(0.27); // scaling constant of respiration
double qp(0.12); // scaling constant of respiration for producers
double J(8*q); // scaling constant of maximum ingestion rate
double a(0.1); // scaling constant of attack rate <0.5
// species constants
double e[3]{0,0.45,0.85}; // assimilation efficiency of detritus, primary producers and consumers
double CN[3]{0,6.6,5}; // C:N ratio of detritus, primary producers and consumers
double B(0.001); // density dependence [0.001 - 0.1]
double K(10); // // half saturation constant of nutrients uptake
double FR(1); // Hill exponent for the functional response
double Mmax(1); // maximal range of species sizes
double Mmin(-5); // size of the smallest organism

double t = 0.0, tFinal = 10000; // time span of integration
double tRecord = 9000; // time from witch recording begins
double tStep = 1; // time step
double h = 1e-6; // absolute accuracy
double extinctionThreshold (pow(10,-30)); // extinction biomass threshold
double interactionThreshold (0.01); // foraging effort from which a link exists

Community foodWeb(diversity
                  ,nbPrimaryProducer
                  ,nbResource
                  ,dim
                  ,I
                  ,L
                  ,b
                  ,delta
                  ,d
                  ,A
                  ,B
                  ,FR);
foodWeb.setMass(Mmin, Mmax); // create the vector of body mass of organisms sorted by ascendant order
foodWeb.setObject(); // create the vectors of populations and parameters
foodWeb.setSpecies(q, qp, B, e, CN, r, K, J, a, b, P0); // create the primary producers
foodWeb.setForagingEffort(); // create the foraging effort and set the interactions
foodWeb.setResources(R0); // initialize the vector of populations for the resources
foodWeb.Dynamic(t, tFinal, tRecord, tStep, h, extinctionThreshold, interactionThreshold, path);
foodWeb.output(path);

    return 0;
}
