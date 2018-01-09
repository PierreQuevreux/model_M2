#ifndef COMMUNITY_H_INCLUDED
#define COMMUNITY_H_INCLUDED

#include "species.h"
#include "string"
using namespace std;

class Community
{
private :
    //PARAMETERS
    int _diversity; // number of species
    int _nbPrimaryProducer; // number of primary producer
    int _nbResource; // number of resources
    int _dim; // dimension of the system (resources+species)
    double _I; // nutrient input
    double _L; // nutrient leaching rate
    double _b; // prey/predator body mass ratio limit
    double _delta; // fraction of mineralized excreted nutrients
    double _d; // intrinsic mineralization rate of detritus
    double _A; // adaptability
    double _B; // density dependence
    double _FR; // Hill exponent for the functional response
    double _Mmax; // maximal size of organisms
    double _Mmin; // size of the smallest organism
    double _R0; // initial quantity of nutrients
    double _D0; // initial quantity of detritus

    string _output;

    //STRUCTURES
    double *_M; // pointer toward the vector of body mass

    //INTEGRATION STRUCTURES
    int _nba; // number of possible interactions
    double *_pop; // pointer toward the vector of R,D,X and a
    double *_popRecy; // pointer toward the vector of nutrient recycled by each species
    Species **_params; // pointer toward the vector of pointers towards the instanced classes
    double *_TL; // pointer toward the vector of trophic levels
    double *_meanTL; // pointer toward the vector of trophic levels
    double *_tExt; // pointer toward the vector of extinction time of species
    double *_biomass; // sum of biomass for each species
    double *_biomassSQ; // sum of square of biomass for each species
    double *_CV; // coefficient of variation of the biomass of species
    int _nPoint; // number of recorded points counter

    //RECORDED OUT PUT
    double _flux[3]; // primary production, secondary production and recycled nutrients
    double _sumFlux[3]; // total primary production, secondary production and recycled nutrients
    double _sumFluxSQ[3]; // total square primary production, secondary production and recycled nutrients
    double _fluxCV[3]; //coefficient of variation of primary production, secondary production and recycled nutrients
    double _recy[2]; // nutrients recycled by primary and secondary producers
    double _connectance;
    double _TLmax; // maximal TL
    int _diversityFinal; // number of species at the end
    int _nbPrimaryProducerFinal; // number of primary producer at the end

public :
    // CREATION OF THE FOODWEB
    void setMass(double Mmin, double Mmax); // create the vector of body mass of organisms sorted by ascendant order
    void setObject(); // create the vectors of populations and parameters
    void setSpecies(double q, double qp, double B, double e[], double CN[], double r, double K, double J, double f, double b, double P0[]); // create the primary producers
    void setForagingEffort(); // create the foraging effort and set the interactions
    void setResources(double R0[]); // initialize the vector of populations for the resources
    //INTEGRATION
    int Dynamic(double t, double tFinal, double tRecord, double tStep, double h, double extinctionThreshold, double interactionThreshold, string path);
    void output(string path); // write the results in the recording structures

    //constructor
    Community(int diversity
              ,int nbPrimaryProducer
              ,int nbResources
              ,int dim
              ,double I
              ,double L
              ,double b
              ,double delta
              ,double d
              ,double A
              ,double B
              ,double FR);
    // destructor
    ~Community();

    //INTEGRATION
    friend int func(double t, const double y[], double f[], void *params); // system of ODE

};

#endif // COMMUNITY_H_INCLUDED
