#ifndef SPECIES_H_INCLUDED
#define SPECIES_H_INCLUDED

#include <iostream>
#include "string"
#include <cstdlib> // to use calloc
using namespace std;

class Species
{
protected :
    // GENERAL PARAMETERS
    double _M; // body mass
    double _Q; // loss rate (respiration and death)
    double _B; // density dependency
    double _e; // assimilation efficiency by predators
    double _CN; // proportion of nutrients in the biomass
    int _diversity; // number of species
    int _nbResource; // number of resources
    int _nba; // number of possible interactions
    string _type; // type of the species

    // PRIMARY PRODUCER
    double _K; // half saturation constant of nutrients uptake
    double _r; // maximum growth rate

    // CONSUMER
    double _J; // maximal ingestion rate
    double _a; // attack rate
    double _b; // prey/predator body mass ratio limit

    // PREYS
    int _nbPrey; // number of preys
    Species** _prey; // array containing the pointers toward the preys
    int* _IDprey; // array containing the number of the prey
    int* _wIDprey; // array containing the number of the aij in the vector _pop
    double* _h; // array containing the handling time of preys

    // PREDATORS
    int _nbPred; // number of predators
    Species** _pred; // array containing the pointers toward the predators and the associated aij
    int* _IDpred; // array containing the number of the predators
    int* _wIDpred; // array containing the number of the aij in the vector _pop
    int _countPred; // counter to fill in the array of predators' pointers

public :

    // GENERAL FUNCTIONS
    Species* getAddress(); // return the pointer this
    double getMass(); // function returning the body mass
    string getType(); // function returning the type

    // FUNCTIONS RELATIVE TO PREDATORS
    int getNbPred(); // function returning the number of predators
    void changeNbPred(int nbPred); // change the total number of predators
    void setPred(); // create the vector of predators'pointers
    void addPred(int a,int nPred,Species* prey); // add the pointer of the predator to the predators list
    int getwIDpred(int i); // return the index of the aij

    // FUNCTIONS RELATIVE TO PREYS
    int getNbPrey(); // return the number of preys
    void changeNbPrey(int nbPrey); // change the number of prey
    void setPrey(); // create the array containing the pointers of preys
    void addPrey(int j,int nPrey,Species* prey,int nTotPrey); // add the pointer of a prey to the diet
    void setHandlingTime(int nPrey,Species* prey); // function calculating the handling time
    int getIDprey(int i); // return the index of the ith prey
    int getwIDprey(int i); // return the index of the aij

    //constructor
    Species(double M, double q, double B, double e, double CN, int diversity, int nbResource, int nba, string type, double K, double r, double J, double a, double b);
    // destructor
    ~Species();

//INTEGRATION
friend int func(double t, const double y[], double f[], void *params); // system of ODE

};

#endif // SPECIES_H_INCLUDED
