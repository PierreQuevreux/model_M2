#include "community.h"
#include "functions.h"

#include <iostream>
#include <fstream>
#include <sstream> // to use istringstream and convert string to double
#include "string"
#include "math.h"

using namespace std;

int main()
{
string path("data/");
int nSim(2); // number written on the file (starts at zero)
string N; // string version of the int
N = to_string(nSim); // string version of the int

// Load parameters
string file_path(path + "parameters_data_" + N + ".txt"); // file containing the parameters
bool header(true); // is there a header in the parameters file ?
int n_start(0); // number of the first simulation (useful if the loop was interrupted)
int n_simu(0); // number of simulations
n_simu = load_parameters(file_path,1);
int n_params(0); // number of parameters
n_params = load_parameters(file_path,2);

double **params; // table containing the parameters feeding the simulation
params = new double*[n_simu]; // create the table
create_table(params,n_simu,n_params); // create the table
file_path = path + string("parameters_") + N + string(".txt"); // file containing the parameters
load_table(file_path,n_simu,params,header); // read the file containing the parameters
// modified parameters
int diversity(4); // initial number of species
int nbResource(2); // number of resources
int dim(nbResource+diversity); // dimension of the system
int nbFlux(9); // number of recorded flux (PP;SP;Irecy;RecyPP;RecySP;RecyInd;RecyDir)
//double Imin(1); // minimal I for the bifurcation
//double Imax(400); // maximum I for the bifurcation
//double Istep(3); // step of parameter exploration
//int n(static_cast<int>((Imax-Imin)/Istep)+1); // number of steps
//double I[n]; // table containing the values of I
int nPoint[3*n_simu]; // number of time steps per simulation
for(int i=0; i<n_simu; i++){
//    I[i] = Imin + i * Istep;
    nPoint[3*i] = 0; // NC model
    nPoint[3*i+1] = 0; // C model
    nPoint[3*i+2] = 0; // SC model
}
//double I(20); // nutrient input
double L(0.2); // nutrient leaching rate (arbitrary)
//double delta(0.8); // fraction of direct recycling
//double d(0.2); // detritus decomposition rate
// allometric scaling constants
double r(0.87);  // scaling constant of maximal growth rate
double q(0.27); // scaling constant of respiration (metabolic rate)
double qp(0.12); // scaling constant of respiration for primary producers
double J(8); // scaling constant of maximum ingestion rate
double a(0.1); // scaling constant of attack rate <0.5
double b(0.05); // prey-predator body mass ratio
// species constants
double B(0.001); // scaling constant of density dependent mortality [0.001 - 0.1]
double K(10); // half saturation constant of nutrients uptake
double FR(1); // Hill exponent for the functional response
// integration parameters
double t = 0.0, tFinal = 2000; // time span of integration
double tRecord = 1000; // time from which recording begins
double tStep = 0.1;
double h = 1e-6; // absolute accuracy
double extinctionThreshold (pow(10,-30)); // extinction biomass threshold

/////////////////////
// RECORDING FILES //
/////////////////////
// time series
string monfichier(path + "time_series_" + N + ".txt");
ofstream fileTS (monfichier.c_str()); // creation of the file containing species dynamics
// MEAN
monfichier = path + "mean_" + N + ".txt";
ofstream fileMean (monfichier.c_str()); // creation of the file containing biomass and flows mean values
monfichier = path + "CV_" + N + ".txt";
ofstream fileCV (monfichier.c_str()); // creation of the file containing biomass and flows mean values
//int nparams(3); // number of parameters (I, delta, d)
int nvariables(nbResource+4*diversity+6); // number of recorded variables
int ndata(1+n_params+1+nvariables); // number of variables (+1=t)
string parameterNames("simu_ID;I;delta;d;model"); // name of the parameters of the simulation
string parameterValue; // value of the parameters
string model[3]{"NC","C","SC"}; // no cycling, cycling and simulated cycling
// write the header
fileTS << "t;" << parameterNames; // parameters
fileMean << parameterNames; // parameters
fileCV << parameterNames; // parameters
fileTS << ";N;D"; // resources
fileMean << ";N;D"; // resources
fileCV << ";N;D"; // resources
for (int i=0; i<diversity; i++){
    fileTS << ";X" << i+1; // species biomass
    fileMean << ";X" << i+1; // species biomass
    fileCV << ";X" << i+1; // species biomass
}
for (int i=0; i<diversity; i++){
    fileTS << ";GD_X" << i+1; // growth:death ratio
    fileMean << ";GD_X" << i+1; // growth:death ratio
    fileCV << ";GD_X" << i+1; // growth:death ratio
}
for (int i=0; i<diversity; i++){
    fileTS << ";Prod_X" << i+1; // biomass production
    fileMean << ";Prod_X" << i+1; // biomass production
    fileCV << ";Prod_X" << i+1; // biomass production
}
for (int i=0; i<diversity; i++){
    fileTS << ";Recy_X" << i+1; // direct recycled nutrients
    fileMean << ";Recy_X" << i+1; // direct recycled nutrients
    fileCV << ";Recy_X" << i+1; // direct recycled nutrients
}
fileTS << ";IndRecy;SProd;SRrecy;TProd;TDirRecy;TRecy" << endl; // general flows
fileMean << ";IndRecy;SProd;SRrecy;TProd;TDirRecy;TRecy" << endl; // general flows
fileCV << ";IndRecy;SProd;SRrecy;TProd;TDirRecy;TRecy" << endl; // general flows

/////////////////
// SIMULATIONS //
/////////////////
for(int i=n_start; i<n_simu; i++){
    Community FoodWeb(params[i][0],
                      diversity,
                      nbResource,
                      dim,
                      nbFlux,
                      nvariables,
                      params[i][1],
                      L,
                      0, // delta
                      0, // d
                      K,
                      FR);
    FoodWeb.setSystem(r,qp,q,J,B,a,b); // create all the vectors defining the populations
    // First run without nutrient cycling
    parameterValue = to_string(params[i][0]) + ";" +
                         to_string(params[i][1]) + ";" +
                         to_string(params[i][2]) + ";" +
                         to_string(params[i][3]) + ";" +
                         model[0];
    nPoint[3*i] = FoodWeb.Dynamic(t,tFinal,tRecord,tStep,h,extinctionThreshold,
                                  parameterValue,fileTS,fileMean,fileCV);
    // Second run with nutrient cycling
    FoodWeb.resetVariables(params[i][2],params[i][3],false);
    parameterValue = to_string(params[i][0]) + ";" +
                         to_string(params[i][1]) + ";" +
                         to_string(params[i][2]) + ";" +
                         to_string(params[i][3]) + ";" +
                         model[1];
    nPoint[3*i+1] = FoodWeb.Dynamic(t,tFinal,tRecord,tStep,h,extinctionThreshold,
                                    parameterValue,fileTS,fileMean,fileCV);
    // Final run without nutrient cycling but with the enrichment correction
    FoodWeb.resetVariables(0,0,true);
    parameterValue = to_string(params[i][0]) + ";" +
                         to_string(params[i][1]) + ";" +
                         to_string(params[i][2]) + ";" +
                         to_string(params[i][3]) + ";" +
                         model[2];
    nPoint[3*i+2] = FoodWeb.Dynamic(t,tFinal,tRecord,tStep,h,extinctionThreshold,
                                    parameterValue,fileTS,fileMean,fileCV);
    cout << "simulation " << i+1 << "/" << n_simu << " - COMPLETED" << endl;
}
fileTS.close();

/////////////////
// BIFURCATION //
/////////////////
//int nrow(0); // total number of recorded points
//for(int i=0; i<3*n_simu; i++){
//    nrow += nPoint[i];
//}
//cout << nrow << " points for bifurcation" << endl;
//double** TS = new double*[nrow]; // table containing all time series
//create_table(TS,nrow,ndata);
//
//// OPENING OF THE FILE AND STORAGE OF THE INFORMATION IN THE TABLE
//load_table(monfichier,nrow,TS,true);
//
//// FINDING OF EXTREMA AND WRITING IN THE BIFURCATION FILE
//// BIFURCATION
//monfichier = path + "bifurcation_" + N + ".txt";
//ofstream fileBifurcation (monfichier.c_str());
//fileBifurcation << parameterNames << "value;variable" << endl; // value (minimum or maximum), variable ((recorded variable)
//
//// MEAN
//monfichier = path + "mean_" + N + ".txt";
//ofstream fileMean (monfichier.c_str()); // creation of the file containing biomass and flows mean values
//// write the header
//fileMean << parameterNames; // parameters
//fileMean << ";N;D"; // resources
//for (int i=0; i<diversity; i++){
//    fileMean << ";X" << i+1; // species biomass
//}
//for (int i=0; i<diversity; i++){
//    fileMean << ";GD_X" << i+1; // growth:death ratio
//}
//for (int i=0; i<diversity; i++){
//    fileMean << ";Prod_X" << i+1; // biomass production
//}
//for (int i=0; i<diversity; i++){
//    fileMean << ";Recy_X" << i+1; // direct recycled nutrients
//}
//fileMean << ";IndRecy;SProd;SRrecy;TProd;TDirRecy;TRecy" << endl; // general flows
//
//// COEFFICIENT OF VARIATION
//monfichier = path + "CV_" + N + ".txt";
//ofstream fileCV (monfichier.c_str()); // creation of the file containing biomass and flows mean values
//fileCV << parameterNames; // parameters
//fileCV << ";N;D"; // resources
//for (int i=0; i<diversity; i++){
//    fileCV << ";X" << i+1; // species biomass
//}
//for (int i=0; i<diversity; i++){
//    fileCV << ";GD_X" << i+1; // growth:death ratio
//}
//for (int i=0; i<diversity; i++){
//    fileCV << ";Prod_X" << i+1; // biomass production
//}
//for (int i=0; i<diversity; i++){
//    fileCV << ";Recy_X" << i+1; // direct recycled nutrients
//}
//fileCV << ";IndRecy;SProd;SRrecy;TProd;TDirRecy;TRecy" << endl; // general flows
//
//// Perform the bifurcation
//int nstart(0); // starting point in the time series table
//for (int i=0; i<3*n_simu; i++){
//    bifurcation(fileBifurcation,variables,nPoint[i],nstart,n_params,nvariables,TS);
//    stat(fileMean,fileCV,variables,nPoint[i],nstart,n_params,nvariables,TS);
//    nstart += nPoint[i];
//    cout << "bifurcation " << i+1 << "/" << 3*n_simu << " - COMPLETED" << endl;//}
//
// Delete the table of parameters
delete_table(params,n_simu);
//delete_table(TS,nrow);

    return 0;

}
