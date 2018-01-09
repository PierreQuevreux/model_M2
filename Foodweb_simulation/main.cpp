#include "community.h"
#include "species.h"
#include "functions.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include "string"
#include "math.h"

using namespace std;

int main()
{
string path("results/"); // path to the folder where results are saved
stringstream out; // to convert int into string
string N; // string version of the int

// VARIABLES
int nSim(0); // number written on the file (starts at zero)
int nI(5); // number of I values tested
//double Ival[50]{1, 3, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400}; // variable testes
double Ival[5]{1, 25, 110, 210, 310}; // variable testes
//double Ival[5]{3, 27.5, 120, 220, 320}; // variable testes
//double Ival[5]{5, 30, 130, 230, 330}; // variable testes
//double Ival[5]{7.5, 40, 140, 240, 340}; // variable testes
//double Ival[5]{10, 50, 150, 250, 350}; // variable testes
//double Ival[5]{12.5, 60, 160, 260, 360}; // variable testes
//double Ival[5]{15, 70, 170, 270, 370}; // variable testes
//double Ival[5]{17.5, 80, 180, 280, 380}; // variable testes
//double Ival[5]{20, 90, 190, 290, 390}; // variable testes
//double Ival[5]{22.5, 100, 200, 300, 400}; // variable testes
int replicate(50); // number of replicate for each variable set

// PARAMETERS
// Ecosystem
int diversity(50); // initial number of species
int nbPrimaryProducer(5); // initial number of primary producers
int nbResource(2); // number of resources
int dim(nbResource+diversity); // dimension of the system
//double I(10); // nutrient input
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
double a(0.1); // scaling constant of attack rate
// species constants
double e[3]{0,0.45,0.85}; // assimilation efficiency of detritus, primary producers and consumers
double CN[3]{0,6.6,5}; // C:N ratio of detritus, primary producers and consumers
double B(0.001); // density dependence (0.00001 max)
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

// FILE WITH GENERAL DATA
out << nSim;
N = out.str();

string monfichier(path + "data" + N + ".txt");
ofstream filedata (monfichier.c_str()); // creation of the file of populations
filedata << "nbSimu" <<
            ";I" <<
            ";delta" <<
            ";d" <<
            ";NbSpeciesInit" <<
            ";NbPPInit" <<
            ";NbSpeciesFinal" <<
            ";NbPPFinal" <<
            ";PP" <<
            ";SP" <<
            ";Irecy" <<
            ";PPcv" <<
            ";PScv" <<
            ";Irecycv" <<
            ";TLmax" <<
            ";C" <<
            ";CV" <<
            ";RecyPP" <<
            ";RecySP" << endl; // write the header


// FILE WITH SPECIES AVERAGE BIOMASSES
monfichier = path + "species" + N + ".txt";
ofstream filespecies (monfichier.c_str()); // creation of the file of populations
filespecies << "nbSimu;I;delta;d;x1" ; // write the header
for (int i=2; i<=diversity; i++){
    filespecies << ";x" << i ;
}
filespecies << endl;

// FILE WITH SPECIES AVERAGE COEFFICIENT OF VARIATION
monfichier = path + "CV" + N + ".txt";
ofstream fileCV (monfichier.c_str()); // creation of the file of populations
fileCV << "nbSimu;I;delta;d;x1" ; // write the header
for (int i=2; i<=diversity; i++){
    fileCV << ";x" << i ;
}
fileCV << endl;

// FILE WITH SPECIES AVERAGE TROPHIC LEVEL
monfichier = path + "TL" + N + ".txt";
ofstream fileTL (monfichier.c_str()); // creation of the file of populations
fileTL << "nbSimu;I;delta;d;x1" ; // write the header
for (int i=2; i<=diversity; i++){
    fileTL << ";x" << i ;
}
fileTL << endl;

int ndata(19); // number of recorded results
double data[ndata]; // array containing the results
double species[diversity]; // array with the average biomass of species
double CV[diversity]; // array with the average coefficient of variation of species
double TL[diversity]; // array with the average trophic level of species
double nbSimu(0); // counter of the number of simulations

for (int i=0; i<nI; i++){

    for (int j=0; j<replicate; j++){
        // Variables
        double I = Ival[i]; // density dependence (0.00001 max)

        // Simulation
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
        foodWeb.setResources(R0); // initialise the vector of populations for the resources

    /////////////////////////////////
    // FIRST RUN WITHOUT RECYCLING //
    /////////////////////////////////
        foodWeb.resetVariables(I, 0, 0); // remove nutrient recycling
        foodWeb.Dynamic(t, tFinal, tRecord, tStep, h, extinctionThreshold, interactionThreshold, false);
        foodWeb.output(data, species, CV, TL);

        nbSimu++;
        data[0] = nSim * nI * replicate * 3 + nbSimu;
        data[1] = I;
        data[2] = 0; // delta
        data[3] = 0; // d
        // Writing in the data file
        for (int k=0; k<ndata-1; k++){
            filedata << data[k] << ";";
        }
        filedata << data[ndata-1] << endl;

        cout << nbSimu << "/" << nI * replicate * 3 << endl;

        // Writing in the species file
        filespecies << data[0] << ";"; // number of the simulation
        filespecies << data[1] << ";"; // I
        filespecies << data[2] << ";"; // delta
        filespecies << data[3] << ";"; // d
        for (int k=0; k<diversity-1; k++){ // species' biomass
            filespecies << species[k] << ";";
        }
        filespecies << species[diversity-1] << endl;

        // Writing in the CV file
        fileCV << data[0] << ";"; // number of the simulation
        fileCV << data[1] << ";"; // I
        fileCV << data[2] << ";"; // delta
        fileCV << data[3] << ";"; // d
        for (int k=0; k<diversity-1; k++){ // species' CV
            fileCV << CV[k] << ";";
        }
        fileCV << CV[diversity-1] << endl;

        // Writing in the TL file
        fileTL << data[0] << ";"; // number of the simulation
        fileTL << data[1] << ";"; // I
        fileTL << data[2] << ";"; // delta
        fileTL << data[3] << ";"; // d
        for (int k=0; k<diversity-1; k++){ // species' TL
            fileTL << TL[k] << ";";
        }
        fileTL << TL[diversity-1] << endl;

    /////////////////////////////////////////////////
    // SECOND RUN OF THE SIMULATION WITH RECYCLING //
    /////////////////////////////////////////////////
        foodWeb.resetVariables(I, delta, d); // remove nutrient recycling and correct the nutrient input
        foodWeb.Dynamic(t, tFinal, tRecord, tStep, h, extinctionThreshold, interactionThreshold, true);
        foodWeb.output(data, species, CV, TL);

        nbSimu++;
        data[0] = nSim * nI * replicate * 3 + nbSimu;
        data[1] = I;
        data[2] = delta;
        data[3] = d;
        // Writing in the data file
        for (int k=0; k<ndata-1; k++){
            filedata << data[k] << ";";
        }
        filedata << data[ndata-1] << endl;

        cout << nbSimu << "/" << nI * replicate * 3 << endl; // number of executed simulation to the total number of simulations to run

        // Writing in the species file
        filespecies << data[0] << ";"; // number of the simulation
        filespecies << data[1] << ";"; // I
        filespecies << data[2] << ";"; // delta
        filespecies << data[3] << ";"; // d
        for (int k=0; k<diversity-1; k++){ // species' biomass
            filespecies << species[k] << ";";
        }
        filespecies << species[diversity-1] << endl;

        // Writing in the CV file
        fileCV << data[0] << ";"; // number of the simulation
        fileCV << data[1] << ";"; // I
        fileCV << data[2] << ";"; // delta
        fileCV << data[3] << ";"; // d
        for (int k=0; k<diversity-1; k++){ // species' CV
            fileCV << CV[k] << ";";
        }
        fileCV << CV[diversity-1] << endl;

        // Writing in the TL file
        fileTL << data[0] << ";"; // number of the simulation
        fileTL << data[1] << ";"; // I
        fileTL << data[2] << ";"; // delta
        fileTL << data[3] << ";"; // d
        for (int k=0; k<diversity-1; k++){ // species' TL
            fileTL << TL[k] << ";";
        }
        fileTL << TL[diversity-1] << endl;


    ///////////////////////////////////////////////////////////////////
    // THIRD RUN WITHOUT RECYCLING BUT WITH CORRECTED NUTRIENT INPUT //
    ///////////////////////////////////////////////////////////////////
        foodWeb.resetVariables(I+data[10], 0, 0); // remove nutrient recycling and correct the nutrient input
        foodWeb.Dynamic(tRecord, tFinal, tRecord, tStep, h, extinctionThreshold, interactionThreshold, false);
        foodWeb.output(data, species, CV, TL);

        nbSimu++;
        data[0] = nSim * nI * replicate * 3 + nbSimu;
        data[1] = I;
        data[2] = -1; // delta
        data[3] = -1; // d
        // Writing in the data file
        for (int k=0; k<ndata-1; k++){
            filedata << data[k] << ";";
        }
        filedata << data[ndata-1] << endl;

        cout << nbSimu << "/" << nI * replicate * 3 << endl;

        // Writing in the species file
        filespecies << data[0] << ";"; // number of the simulation
        filespecies << data[1] << ";"; // I
        filespecies << data[2] << ";"; // delta
        filespecies << data[3] << ";"; // d
        for (int k=0; k<diversity-1; k++){ // species' biomass
            filespecies << species[k] << ";";
        }
        filespecies << species[diversity-1] << endl;

        // Writing in the species file
        fileCV << data[0] << ";"; // number of the simulation
        fileCV << data[1] << ";"; // I
        fileCV << data[2] << ";"; // delta
        fileCV << data[3] << ";"; // d
        for (int k=0; k<diversity-1; k++){ // species' CV
            fileCV << CV[k] << ";";
        }
        fileCV << CV[diversity-1] << endl;

        // Writing in the species file
        fileTL << data[0] << ";"; // number of the simulation
        fileTL << data[1] << ";"; // I
        fileTL << data[2] << ";"; // delta
        fileTL << data[3] << ";"; // d
        for (int k=0; k<diversity-1; k++){ // species' TL
            fileTL << TL[k] << ";";
        }
        fileTL << TL[diversity-1] << endl;

    }
}
    return 0;
}
