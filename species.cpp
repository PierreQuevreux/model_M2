#include "species.h"
#include "functions.h"

#include <iostream>
#include "string"
#include "math.h"
#include <cstdlib> // to use calloc
using namespace std;

// GENERAL FUNCTIONS
Species* Species::getAddress(){
    return (this);
}
double Species::getMass(){
    return(_M);
}
string Species::getType(){
    return(_type);
}

// FUNCTIONS RELATIVE TO PREDATORS
int Species::getNbPred(){
    return(_nbPred);
}
void Species::changeNbPred(int nbPred){
    _nbPred = nbPred;
}
void Species::setPred(){
    _pred = new Species*[_nbPred];
        for (int i=0; i<_nbPred; i++){
            _pred[i] = NULL;
        }
    _IDpred = new int [_nbPred]; // create the vector of predators' numbers
    _wIDpred = new int [_nbPred]; // create the vector of aij' numbers
}
void Species::addPred(int a,int nPred,Species* pred){
    _pred[_countPred] = pred; // pointer toward the predator
    _IDpred[_countPred] = nPred; // position of the predator
    _wIDpred[_countPred] = _nbResource + _diversity + a; //position of the associated aij
    _countPred++; // counter of predators
}
int Species::getwIDpred(int i){
    return (_wIDpred[i]);
}

// FUNCTIONS RELATIVE TO PREYS
int Species::getNbPrey(){
    return(_nbPrey);
}
void Species::changeNbPrey(int nbPrey){
    _nbPrey=nbPrey;
}
void Species::setPrey(){
    _prey = new Species*[_nbPrey];
        for (int i=0; i<_nbPrey; i++){
            _prey[i] = NULL;
        }
    _h = new double [_nbPrey]; // create the vector of handling times
    _IDprey = new int [_nbPrey]; // create the vector of preys' numbers
    _wIDprey = new int [_nbPrey]; // create the vector of aij' numbers
}
void Species::addPrey(int j,int nPrey,Species* prey,int nTotPrey){
    _prey[nPrey] = prey; // pointer toward the prey
    _IDprey[nPrey] = _nbResource + j; // position of the prey in the vector
    _wIDprey[nPrey] = _nbResource + _diversity + nTotPrey; // nTotPrey : total of preys for all predators, position of the associated aij
}
void Species::setHandlingTime(int nPrey,Species* prey){
    //_h[nPrey]=_b/(2*_J*(_b - (prey->getMass()/_M)));
    _h[nPrey] = _b*_b/(6*_J*(_b - (prey->getMass()/_M)))*(_M/prey->getMass());
}
int Species::getIDprey(int i){
    return (_IDprey[i]);
}
int Species::getwIDprey(int i){
    return (_wIDprey[i]);
}

// constructor
Species::Species(double M
                 ,double q
                 ,double B
                 ,double e
                 ,double CN
                 ,int diversity
                 ,int nbResource
                 ,int nba
                 ,string type
                 ,double K
                 ,double r
                 ,double J
                 ,double a
                 ,double b)
                 :_M(M)
                 ,_e(e)
                 ,_CN(CN)
                 ,_diversity(diversity)
                 ,_nbResource(nbResource)
                 ,_nba(nba)
                 ,_type(type)
                 ,_b(b)
                 ,_nbPrey(0)
                 ,_nbPred(0)
                 ,_countPred(0){
                     if(type == "Plant"){
                         _K = K;
                         _Q = setAllometric(q,M);
                         _B = setAllometric(B,M);
                         _J = 0;
                         _a = 0;
                         _r = setAllometric(r,M);
                     }
                     else{
                        _K = 0;
                        _Q = setAllometric(q,M);
                        _B = setAllometric(B,M);
                        _J = setAllometric(J,M);
                        _a = setAllometric(a,M);
                        _r = 0;
                     }
                 }
// destructor
Species::~Species(){
    if (_nbPred != 0){
        delete[] _pred;
        delete[] _IDpred;
        delete[] _wIDpred;
    }
    if (_nbPrey != 0){
        delete[] _prey;
        delete[] _IDprey;
        delete[] _wIDprey;
        delete[] _h;
    }
}
