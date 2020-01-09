#ifndef COMMUNITY_H_INCLUDED
#define COMMUNITY_H_INCLUDED

#include "string"
using namespace std;

class Community
{
private :
    int _simu_ID; // number of the simulation
    // PARAMETERS
    int _diversity; // initial number of species
    int _nbResource; // number of resources
    int _dim; // dimension of the system
    int _nbFlux; // number of recorded flux (PP;SP;Irecy;RecyPP;RecySP;RecyInd;RecyDir)
    int _nvariables; // number of recorded variables
    double _I; // nutrient input
    double _L; // nutrient leaching rate (arbitrary)
    double _delta; // fraction of mineralised excreted nutrients
    double _d; // intrinsic mineralisation rate of detritus
    // parameters of the ecosystem
    double _g; // growth rate of the primary producer
    double _K; // primary producers half saturation
    double _FR; // Hill exponent for the functional response
    double* _M; // body masses of species
    double* _x; // metabolic rate
    double* _B; // density dependent mortality rate
    double* _a; // attack rate
    double* _h; // handling time
    double* _e; // assimilation efficiency
    double* _CN; // C:N ratio
    double* _pop; // initial biomasses
    double* _flux; // array with the primary production, herbivore production, carnivore production, primary prod. direct recycling, herbivore direct recycling and carnivore direct recycling and indirect recycling
    double* _GD; // growth/death ratio of species
    double* _mean; // mean of recorded variables
    double* _CV; // CV of recorded variables
    double _TR; // total recycling
    double _meanTR; // average total recycling
    double _Irecy; // correction of the enrichment effect
    bool _corr; // simulation with a correction ?
    // tables for mean and CV

public :
    // CREATION OF THE FOODWEB
    void setSystem(double r, double qp, double q, double J, double B, double a, double b);
    void resetVariables(double delta, double d, bool corr);
    //INTEGRATION
    int Dynamic(double t, double tFinal, double tRecord, double tStep, double h, double extinctionThreshold,
                string parameterValue, ofstream &fileTimeSeries, ofstream &fileMean, ofstream &fileCV);

    //constructor
    Community(
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
        double FR);
    // destructor
    ~Community();

    //INTEGRATION
    friend int func(double t, const double y[], double f[], void *params); // system of ODE
    friend int funcBis(double t, const double y[], double f[], void *params); // system of ODE
};

#endif // COMMUNITY_H_INCLUDED
