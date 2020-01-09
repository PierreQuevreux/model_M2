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
    int _nbFlux; // number of ecosystem flows
    int _nbTotB; // number of recorded aggregated biomasses
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
    Species **_params; // pointer toward the vector of pointers towards the instanced classes
    double *_TL; // pointer toward the vector of trophic levels
    double *_meanTL; // pointer toward the vector of average trophic levels
    double *_tExt; // pointer toward the vector of extinction time of species
    double *_meanBiomass; // sum of biomass for each species
    double *_biomassSQ; // sum of square of biomass for each species
    double *_biomassCV; // coefficient of variation of the biomass of species
    double *_recy; // nutrients recycled by each species
    double *_meanRecy; // average amount of nutrients recycled by each species
    double *_recySQ; // sum of square of nutrients recycled by each species
    double *_recyCV; // coefficient of variation of nutrients recycled by each species
    double *_detritus; // detritus produced by each species
    double *_meanDetritus; // average amount of detritus produced by each species
    double *_detritusSQ; // sum of square of detritus produced by each species
    double *_detritusCV; // coefficient of variation of detritus produced by each species
    double *_flux; // recycled nutrients, nutrients recycled by primary, secondary producers, total direct recycling, indirect recycling
    double *_sumFlux; // total primary production, secondary production, recycled nutrients, nutrients recycled by primary and secondary producers and indirect recycling
    double *_sumFluxSQ; // total square primary production, secondary production, recycled nutrients and nutrients recycled by primary and secondary producers
    double *_fluxCV; //coefficient of variation of primary production, secondary production, recycled nutrients and nutrients recycled by primary and secondary producers
    double *_totalBiomass; // total biomass of all species, PP and SP
    double *_meanTotalBiomass; // sum of total biomass of all species, PP and SP
    double *_totalBiomassSQ; // sum of square of total biomass of all species, PP and SP
    double *_totalBiomassCV; // coefficient of variation of the total biomass of all species, PP and SP
    double *_prod; // biomass production of all species, PP and SP
    double *_meanProd; // sum of biomass production of all species, PP and SP
    double *_prodSQ; // sum of square of biomass production of all species, PP and SP
    double *_prodCV; // coefficient of variation of the biomass production of all species, PP and SP
    double *_mort; // total density dependent mortality of all species, PP and SP
    double *_meanMort; // sum of total density dependent mortality of all species, PP and SP
    double *_mortSQ; // sum of square of total density dependent mortality of all species, PP and SP
    double *_mortCV; // coefficient of variation of the total density dependent mortality of all species, PP and SP
    int _nPoint; // number of recorded points counter

    //RECORDED OUT PUT
    double _connectance;
    double _TLmax; // maximal TL
    int _diversityFinal; // number of species at the end
    int _nbPrimaryProducerFinal; // number of primary producer at the end

public :
    // CREATION OF THE FOODWEB
    void setMass(double Mmin, double Mmax, int seed); // create the vector of body mass of organisms sorted by ascendant order
    void setObject(); // create the vectors of populations and parameters
    void setSpecies(double q, double qp, double B, double e[], double CN[], double r, double K, double J, double f, double b, double **pop); // create the primary producers
    void setForagingEffort(); // create the foraging effort and set the interactions
    void setResources(double **pop); // initialize the vector of populations for the resources
    //INTEGRATION
    void resetVariables(double I, double delta, double d); // change some variables
    int Dynamic(double t, double tFinal, double tRecord, double tStep, double h, double extinctionThreshold, double interactionThreshold, bool recyRec,
                string parameterValue, bool time_series_record, ofstream& fileTimeSeriesData, ofstream& fileTimeSeriesBiomass, ofstream& fileTimeSeriesRecy, ofstream& fileTimeSeriesDetritus);
    void output(double **data, double M[], double **biomass, double **biomassCV, double **recy, double **recyCV, double **detritus, double **detritusCV, double **TL, double **tExt, int run); // write the results in the recording structures

    //constructor
    Community(int diversity
              ,int nbPrimaryProducer
              ,int nbResources
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
              ,double FR);
    // destructor
    ~Community();

    //INTEGRATION
    friend int func(double t, const double y[], double f[], void *params); // system of ODE

};

#endif // COMMUNITY_H_INCLUDED
