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

// modified parameters
double Imin(1); // minimal I for the bifurcation
double Imax(40); // maximum I for the bifurcation
double Istep(0.1); // step of the exploration of I
int n(static_cast<int>((Imax-Imin)/Istep)+1); // number of steps
double I[n]; // carrying capacity
int nPoint[n]; // number of simulations
for(int i=0; i<n; i++){
    I[i] = Imin + i * Istep;
    nPoint[i] = 0;
}
//double I(20); // nutrient input
double L(0.2); // nutrient leaching rate (arbitrary)
double delta(0.2); // fraction of mineralized excreted nutrients 0.2
double d(0.2); // intrinsic mineralization rate of detritus 0.2
// allometric scaling constants
double r(0.87);  // scaling constant of maximal growth rate
double q(0.27); // scaling constant of respiration
double qp(0.12); // scaling constant of respiration for producers
double a(0.07); // scaling constant of attack rate <0.5
double b(0.05); // prey-predator body mass ratio
// species constants
//double B(0.001); // density dependence [0.001 - 0.1]
double B(0); // density dependence [0.001 - 0.1]
double K(10); // half saturation constant of nutrients uptake
double FR(1); // Hill exponent for the functional response
// integration parameters
double t = 0.0, tFinal = 1000; // time span of integration
double tRecord = 900; // time from witch recording begins
double tStep = 1;
double h = 1e-6; // absolute accuracy
double extinctionThreshold (pow(10,-30)); // extinction biomass threshold

/////////////////////
// RECORDING FILES //
/////////////////////
string monfichier(path + "chronic.txt");
ofstream fileChronic (monfichier.c_str()); // creation of the file of populations
string variables[19]={"t","I","delta","d","N","D","P","H","C","PP","HP","CP","PR","HR","CR","SP","SR","TP","TR"};
fileChronic << "t;I;delta;d;N;D;P;H;C;PP;HP;CP;PR;HR;CR;SP;SR;TP;TR" << endl; // write the header

/////////////////
// SIMULATIONS //
/////////////////
for(int i=0; i<n; i++){
    Community FoodWeb(
            I[i]
            ,L
            ,delta
            ,d
            ,r
            ,qp
            ,q
            ,B
            ,K
            ,FR);
    FoodWeb.setSystem(qp,q,a,b); // create the vectors of masses
    nPoint[i] = FoodWeb.Dynamic(t,tFinal,tRecord,tStep,h,extinctionThreshold,fileChronic);
    //cout << i+1 << endl;
}
fileChronic.close();


/////////////////
// BIFURCATION //
/////////////////
int N(0); // total number of recorded points
for(int i=0; i<n; i++){
    N += nPoint[i];
}
cout << N << endl;
double** chronic = new double*[N]; // table with the data (lines)
for (int i=0; i<N; i++){
    chronic[i] = new double[19]; // columns
}

// OPENING OF THE FILE AND STORAGE OF THE INFORMATION IN THE TABLE
ifstream file(monfichier.c_str());
if (file){}
else {cout << "ERROR: unable to open file at path" << endl;}
string line;
string car;
getline(file, line); // skip the header
int k(0); // counter for the table of names
string var[19]; // string array containing the string form of the parameters
istringstream iss; // variable for the conversion from string to double
for(int i=0; i<N; i++){
    getline(file, line);
    k=0;
    for (int j=0; j<19; j++){
        var[j].clear(); // clear the string array
    }
    for (unsigned int j=0; j<line.size(); j++){
        car = line[j];
        if (car!=";"){
            var[k] += line[j]; // write the data letter by letter
        }
        else {k++;}
    }
    for (int j=0; j<19; j++){
        iss.str( var[j].c_str() ); // K
        iss >> chronic[i][j];
        iss.clear();
    }
}

// FINDING OF EXTREMA AND WRITING IN THE BIFURCATION FILE
monfichier = path + "bifurcation.txt";
ofstream fileBifurcation (monfichier.c_str()); // creation of the file of populations
fileBifurcation << "I;delta;d;value;variable" << endl; // write the header
k=0;
for(int i=1; i<nPoint[0]-1; i++){
    for(int j=4; j<19; j++){
        if( (chronic[i][j]<=chronic[i-1][j] && chronic[i][j]<=chronic[i+1][j]) || (chronic[i][j]>=chronic[i-1][j] && chronic[i][j]>=chronic[i+1][j]) ){
            fileBifurcation << chronic[i][1] << ";" << chronic[i][2] << ";" << chronic[i][3] << ";" << chronic[i][j] << ";" << variables[j] << endl;
        }
    }
}
for(int z=1; z<n; z++){
    k += nPoint[z-1];
    for(int i=k+1; i<k+nPoint[z]-1; i++){
        for(int j=4; j<19; j++){
            if( (chronic[i][j]<=chronic[i-1][j] && chronic[i][j]<=chronic[i+1][j]) || (chronic[i][j]>=chronic[i-1][j] && chronic[i][j]>=chronic[i+1][j]) ){
                fileBifurcation << chronic[i][1] << ";" << chronic[i][2] << ";" << chronic[i][3] << ";" << chronic[i][j] << ";" << variables[j] << endl;
            }
        }
    }
}

for (int i=0; i<N; i++){
    delete chronic[i];
}
delete[] chronic;

    return 0;

}
