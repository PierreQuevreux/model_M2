#ifndef COMMUNITY_H_INCLUDED
#define COMMUNITY_H_INCLUDED

#include "string"
using namespace std;

class Community
{
private :
    // PARAMETERS
    double _I; // nutrient input
    double _L; // nutrient leaching rate (arbitrary)
    double _delta; // fraction of mineralized excreted nutrients
    double _d; // intrinsic mineralization rate of detritus
    // parameters of the ecosystem
    double _r; // growth rate of the primary producer
    double _qp; // scaling constant of respiration for producers
    double _q; // scaling constant of respiration
    double _B; // density dependence [0.001 - 0.1]
    double _K; // primary producers half saturation
    double _FR; // Hill exponent for the functional response
    double _M[5]; // body masses of species
    double _x[5]; // initial metabolic rate
    double _a[5]; // initial metabolic rate
    double _h[5]; // initial metabolic rate
    double _e[5]; // assimilation efficiency
    double _CN[5]; // C:N ratio
    double _pop[5]; // initial biomasses
    double _flux[6]; // array with the value of PP and SP at each step

public :
    // CREATION OF THE FOODWEB
    void setSystem(double qp, double q, double a, double b);
    //INTEGRATION
    int Dynamic(double t, double tFinal, double tRecord, double tStep, double h, double extinctionThreshold, ofstream &fileChronic);

    //constructor
    Community(
        double I,
        double L,
        double delta,
        double d,
        double r,
        double qp,
        double q,
        double B,
        double K,
        double FR);
    // destructor
    ~Community();

    //INTEGRATION
    friend int func(double t, const double y[], double f[], void *params); // system of ODE

};

#endif // COMMUNITY_H_INCLUDED
