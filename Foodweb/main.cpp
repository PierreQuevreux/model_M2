#include "community.h"
#include "species.h"
#include "functions.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include "string"
#include "math.h"

#include <time.h> // to measure the execution time

#include <omp.h> // OpenMP library for parallel calculation
// add the option -fopenmp to the compiler

using namespace std;

int main()
{
    clock_t tStart = clock(); // to measure the execution time

    int nSim(0); // number written on the file (starts at zero)
    stringstream out; // to convert int into string
    string N; // string version of the int
    out << nSim;
    N = out.str();
    bool time_series_record = true; // record or not the time series

    ////////////////////////////////
    // PARALLELISATION PARAMETERS //
    ////////////////////////////////
    string path("results/"); // path to the folder where results are saved
    string file_path(path + "parameters_data_" + N + ".txt"); // file containing the parameters
    bool header(true); // is there a header in the parameters file ?
    int n_simu(0); // number of simulations
    n_simu = load_parameters_data(file_path,1);
    int n_params(0); // number of parameters
    n_params = load_parameters_data(file_path,2);
    int n_threads(0); // number of threads
    n_threads = load_parameters_data(file_path,3);
    omp_set_num_threads(n_threads); // number of threads used

    int n_start(0); // number of the first simulation (useful if the loop was interrupted)
    double **params; // table containing the parameters feeding the simulation
    params = new double*[n_simu]; // create the table
    create_table(params,n_simu,n_params); // create the table
    file_path = path + string("parameters_") + N + string(".txt"); // file containing the parameters
    load_parameters(file_path,n_simu,params,header); // read the file containing the parameters

    /////////////////////////////
    // PROGRAM FOR SIMULATIONS //
    /////////////////////////////

    // PARAMETERS
    // Ecosystem
    int diversity(50); // initial number of species
    int nbPrimaryProducer(5); // initial number of primary producers
    int nbResource(2); // number of resources
    int dim(nbResource+diversity); // dimension of the system
    int nbFlux(5); // number of recorded flux (Irecy;RecyPP;RecySP;RecyDir;RecyInd)
    int nbTotB(3); // number of recorded aggregated biomasses
    //double I(10); // nutrient input
    double L(0.2); // nutrient leaching rate (arbitrary)
    //double delta(0.2); // fraction of mineralized excreted nutrients
    //double d(0.2); // intrinsic mineralization rate of detritus
    // Initialisation
    //double R0[2]{10,10}; // initial density of nutrients and detritus
    //double P0[2]{10,5}; // initial density of primary producers and predators
    double **pop; // table containing the parameters feeding the simulation
    pop = new double*[dim]; // create the table
    create_table(pop,1,dim); // create the table
    file_path = path + string("initialisation_") + N + string(".txt"); // file containing the parameters
    load_parameters(file_path,1,pop,false); // read the file containing the parameters
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
    double tTrans = 1000; // time to avoid the transient dynamics in the third run =300
    double tStep = 1; // time step
    double h = 1e-6; // absolute accuracy
    double extinctionThreshold (pow(10,-30)); // extinction biomass threshold
    double interactionThreshold (0.01); // foraging effort from which a link exists

    // FILE WITH GENERAL DATA
    int ndata(34); // number of recorded results
    string monfichier(path + "data_" + N + ".txt");
    ofstream fileData (monfichier.c_str()); // creation of the file of populations
    string parameterNames("simu_ID;seed;I;delta;d;model"); // name of the parameters of the simulation
    string model[3]{"NC","C","SC"}; // no cycling, cycling and simulated cycling

    fileData << "NbSpeciesInit;NbPPInit;NbSpeciesFinal;NbPPFinal;TLmax;connectance"; // 6
    fileData << ";totalBiomass;PPbiomass;SPbiomass;totalProd;PPprod;SPprod"; // 6
    fileData << ";Irecy;RecyPP;RecySP;RecyDir;RecyInd;totalMort;PPmort;SPmort"; // 8
    fileData << ";totalBiomassCV;PPbiomassCV;SPbiomassCV;totalProdCV;PPprodCV;SPprodCV"; // 6
    fileData << ";IrecyCV;RecyPPcv;RecySPcv;RecyDirCV;RecyIndCV;totalMortCV;PPmortCV;SPmortCV"; // 8
    fileData << ";" << parameterNames << endl;

    // FILE WITH BODY MASS OF ORGANISM
    monfichier = path + "bodymass_" + N + ".txt";
    ofstream fileMass (monfichier.c_str()); // creation of the file of populations
    fileMass << "nbSimu"; // write the header
    for (int i=1; i<=diversity; i++){
        fileMass << ";x" << i ;
    }
    fileMass << endl;

    // FILE WITH SPECIES AVERAGE BIOMASSES
    monfichier = path + "biomass_" + N + ".txt";
    ofstream fileBiomass (monfichier.c_str()); // creation of the file of populations
    fileBiomass << parameterNames << ";N;D"; // write the header
    for (int i=1; i<=diversity; i++){
        fileBiomass << ";x" << i ;
    }
    fileBiomass << endl;

    // FILE WITH SPECIES AVERAGE COEFFICIENT OF VARIATION
    monfichier = path + "biomassCV_" + N + ".txt";
    ofstream fileBiomassCV (monfichier.c_str()); // creation of the file of populations
    fileBiomassCV << parameterNames << ";N;D"; // write the header
    for (int i=1; i<=diversity; i++){
        fileBiomassCV << ";x" << i ;
    }
    fileBiomassCV << endl;

    // FILE WITH NUTRIENT RECYCLED BY EACH SPECIES
    monfichier = path + "recy_" + N + ".txt";
    ofstream fileRecy (monfichier.c_str()); // creation of the file of populations
    fileRecy << parameterNames; // write the header
    for (int i=1; i<=diversity; i++){
        fileRecy << ";x" << i ;
    }
    fileRecy << endl;

    // FILE WITH THE COEFFICIENT OF VARIATION OF NUTRIENT RECYCLED BY EACH SPECIES
    monfichier = path + "recyCV_" + N + ".txt";
    ofstream fileRecyCV (monfichier.c_str()); // creation of the file of populations
    fileRecyCV << parameterNames; // write the header
    for (int i=1; i<=diversity; i++){
        fileRecyCV << ";x" << i ;
    }
    fileRecyCV << endl;

    // FILE WITH DETRITUS PRODUCED BY EACH SPECIES
    monfichier = path + "detritus_" + N + ".txt";
    ofstream fileDetritus (monfichier.c_str()); // creation of the file of populations
    fileDetritus << parameterNames; // write the header
    for (int i=1; i<=diversity; i++){
        fileDetritus << ";x" << i ;
    }
    fileDetritus << endl;

    // FILE WITH THE COEFFICIENT OF VARIATION OF DETRITUS PRODUCED BY EACH SPECIES
    monfichier = path + "detritusCV_" + N + ".txt";
    ofstream fileDetritusCV (monfichier.c_str()); // creation of the file of populations
    fileDetritusCV << parameterNames; // write the header
    for (int i=1; i<=diversity; i++){
        fileDetritusCV << ";x" << i ;
    }
    fileDetritusCV << endl;

    // FILE WITH SPECIES AVERAGE TROPHIC LEVEL
    monfichier = path + "TL_" + N + ".txt";
    ofstream fileTL (monfichier.c_str()); // creation of the file of populations
    fileTL << parameterNames; // write the header
    for (int i=1; i<=diversity; i++){
        fileTL << ";x" << i ;
    }
    fileTL << endl;

    // FILE WITH SPECIES EXTINCTION TIME
    monfichier = path + "tExtinction_" + N + ".txt";
    ofstream filetExt (monfichier.c_str()); // creation of the file of populations
    filetExt << parameterNames; // write the header
    for (int i=1; i<=diversity; i++){
        filetExt << ";x" << i ;
    }
    filetExt << endl;

    // FILE WITH CHRONICS
    monfichier = path + "time_series_data_" + N + ".txt";
    ofstream fileTimeSeriesData (monfichier.c_str()); // creation of the file of populations
    fileTimeSeriesData << parameterNames << ";t;PP;SP;Irecy;RecyPP;RecySP;RecyInd;RecyDir;totalBiomass;PPbiomass;SPbiomass" << endl; // write the header

    monfichier = path + "time_series_biomass_" + N + ".txt";
    ofstream fileTimeSeriesBiomass (monfichier.c_str()); // creation of the file of populations
    fileTimeSeriesBiomass << parameterNames << ";t;N;D"; // write the header
    for (int i=1; i<=diversity; i++){
        fileTimeSeriesBiomass << ";X" << i ;
    }
    fileTimeSeriesBiomass << endl;

    monfichier = path + "time_series_recy_" + N + ".txt";
    ofstream fileTimeSeriesRecy (monfichier.c_str()); // creation of the file of direct recycling
    fileTimeSeriesRecy << parameterNames << ";t"; // write the header
    for (int i=1; i<=diversity; i++){
        fileTimeSeriesRecy << ";X" << i ;
    }
    fileTimeSeriesRecy << endl;

    monfichier = path + "time_series_detritus_" + N + ".txt";
    ofstream fileTimeSeriesDetritus (monfichier.c_str()); // creation of the file of detritus production
    fileTimeSeriesDetritus << parameterNames << ";t"; // write the header
    for (int i=1; i<=diversity; i++){
        fileTimeSeriesDetritus << ";X" << i ;
    }
    fileTimeSeriesDetritus << endl;

    /////////////////////////
    // LOOP TO PARALLELISE //
    /////////////////////////

    #pragma omp parallel for // instruction to parallelise the following for loop
    for(int i=n_start; i<n_simu; i++){
        // Create tables to record results
        string parameterValue; // value of the parameters
        double **data; // array containing the results
        data = new double*[ndata+n_params]; // create the table
        create_table(data,ndata+n_params,3);
        double M[diversity]; // array with the body mass of species
        double **biomass; // array with the average biomass of species
        biomass = new double*[dim]; // create the table
        create_table(biomass,dim,3);
        double **biomassCV; // array with the average coefficient of variation of species
        biomassCV = new double*[dim]; // create the table
        create_table(biomassCV,dim,3);
        double **recy; // array with the nutrient recycled by each species
        recy = new double*[diversity]; // create the table
        create_table(recy,diversity,3);
        double **recyCV; // array with the coefficient of variation of the nutrients recycled by each species
        recyCV = new double*[diversity]; // create the table
        create_table(recyCV,diversity,3);
        double **detritus; // array with the detritus produced by each species
        detritus = new double*[diversity]; // create the table
        create_table(detritus,diversity,3);
        double **detritusCV; // array with the coefficient of variation of the detritus produced by each species
        detritusCV = new double*[diversity]; // create the table
        create_table(detritusCV,diversity,3);
        double **TL; // array with the average trophic level of species
        TL = new double*[diversity]; // create the table
        create_table(TL,diversity,3);
        double **tExt; // array with the average trophic level of species
        tExt = new double*[diversity]; // create the table
        create_table(tExt,diversity,3);

        int run(0); // number of run of the same simulated food web

        // Simulation
        Community foodWeb(diversity
                          ,nbPrimaryProducer
                          ,nbResource
                          ,dim
                          ,nbFlux
                          ,nbTotB
                          ,params[i][2] // I
                          ,L
                          ,b
                          ,params[i][3] //delta
                          ,params[i][4] //d
                          ,A
                          ,B
                          ,FR);
        foodWeb.setMass(Mmin, Mmax, params[i][1]); // create the vector of body mass of organisms sorted by ascendant order
        //foodWeb.setMass(Mmin, Mmax, seeds[i]); // create the vector of body mass of organisms sorted by ascendant order
        foodWeb.setObject(); // create the vectors of populations and parameters
        foodWeb.setSpecies(q, qp, B, e, CN, r, K, J, a, b, pop); // create species
        foodWeb.setForagingEffort(); // create the foraging effort and set the interactions
        foodWeb.setResources(pop); // initialise the vector of populations for the resources

    /////////////////////////////////
    // FIRST RUN WITHOUT RECYCLING //
    /////////////////////////////////
        run = 0;
        parameterValue = to_string(params[i][0]) + ";" +
                         to_string(params[i][1]) + ";" +
                         to_string(params[i][2]) + ";" +
                         to_string(params[i][3]) + ";" +
                         to_string(params[i][4]) + ";" +
                         model[run];
        foodWeb.resetVariables(params[i][2], 0, 0); // remove nutrient recycling
        foodWeb.Dynamic(t, tFinal, tRecord, tStep, h, extinctionThreshold, interactionThreshold, false,
                        parameterValue, time_series_record, fileTimeSeriesData, fileTimeSeriesBiomass, fileTimeSeriesRecy, fileTimeSeriesDetritus);
        foodWeb.output(data, M, biomass, biomassCV, recy, recyCV, detritus, detritusCV, TL, tExt, run);

        data[0+ndata][run] = params[i][0]; // number of the food web
        data[1+ndata][run] = params[i][1]; // seed
        data[2+ndata][run] = params[i][2]; // I
        data[3+ndata][run] = params[i][3]; // delta
        data[4+ndata][run] = params[i][4]; // d

    /////////////////////////////////////////////////
    // SECOND RUN OF THE SIMULATION WITH RECYCLING //
    /////////////////////////////////////////////////
        run = 1;
        parameterValue = to_string(params[i][0]) + ";" +
                         to_string(params[i][1]) + ";" +
                         to_string(params[i][2]) + ";" +
                         to_string(params[i][3]) + ";" +
                         to_string(params[i][4]) + ";" +
                         model[run];
        foodWeb.resetVariables(params[i][2], params[i][3], params[i][4]); // remove nutrient recycling and correct the nutrient input
        foodWeb.Dynamic(t, tFinal, tRecord, tStep, h, extinctionThreshold, interactionThreshold, true,
                        parameterValue, time_series_record, fileTimeSeriesData, fileTimeSeriesBiomass, fileTimeSeriesRecy, fileTimeSeriesDetritus);
        foodWeb.output(data, M, biomass, biomassCV, recy, recyCV, detritus, detritusCV, TL, tExt, run);

        data[0+ndata][run] = params[i][0]; // number of the food web
        data[1+ndata][run] = params[i][1]; // seed
        data[2+ndata][run] = params[i][2]; // I
        data[3+ndata][run] = params[i][3]; // delta
        data[4+ndata][run] = params[i][4]; // d

    ///////////////////////////////////////////////////////////////////
    // THIRD RUN WITHOUT RECYCLING BUT WITH CORRECTED NUTRIENT INPUT //
    ///////////////////////////////////////////////////////////////////
        run = 2;
        parameterValue = to_string(params[i][0]) + ";" +
                         to_string(params[i][1]) + ";" +
                         to_string(params[i][2]) + ";" +
                         to_string(params[i][3]) + ";" +
                         to_string(params[i][4]) + ";" +
                         model[run];
        foodWeb.resetVariables(params[i][2]+data[12][1], 0, 0); // remove nutrient recycling and correct the nutrient input
        foodWeb.Dynamic(tRecord, tFinal+tTrans, tRecord+tTrans, tStep, h, extinctionThreshold, interactionThreshold, false,
                        parameterValue, time_series_record, fileTimeSeriesData, fileTimeSeriesBiomass, fileTimeSeriesRecy, fileTimeSeriesDetritus); // to have 1000 time steps to avoid the new transient dynamic
        foodWeb.output(data, M, biomass, biomassCV, recy, recyCV, detritus, detritusCV, TL, tExt, run);

        data[0+ndata][run] = params[i][0]; // number of the food web
        data[1+ndata][run] = params[i][1]; // seed
        data[2+ndata][run] = params[i][2]; // I
        data[3+ndata][run] = params[i][3]; // delta
        data[4+ndata][run] = params[i][4]; // d

    ///////////////////////////////
    // WRITE IN THE OUTPUT FILES //
    ///////////////////////////////

        #pragma omp critical
        {
            cout << "simulation " << i << " in thread " << omp_get_thread_num() <<  endl;

            // Writing in the body mass file
            fileMass << data[0+ndata][0] << ";"; // number (ID) of the simulation
            for (int k=0; k<diversity-1; k++){
                fileMass << M[k] << ";";
            }
            fileMass << M[diversity-1] << endl;

            for (int j=0; j<3; j++){

                // Writing in the data file
                for (int k=0; k<ndata+n_params; k++){
                    fileData << data[k][j] << ";";
                }
                fileData << model[j] << endl;

                // Writing in the biomass file
                for (int k=0; k<n_params; k++){
                    fileBiomass << data[ndata+k][j] << ";";
                }
                fileBiomass << model[j] << ";";
                for (int k=0; k<dim-1; k++){
                    fileBiomass << biomass[k][j] << ";";
                }
                fileBiomass << biomass[dim-1][j] << endl;

                // Writing in the biomass CV file
                for (int k=0; k<n_params; k++){
                    fileBiomassCV << data[ndata+k][j] << ";";
                }
                fileBiomassCV << model[j] << ";";
                for (int k=0; k<dim-1; k++){
                    fileBiomassCV << biomassCV[k][j] << ";";
                }
                fileBiomassCV << biomassCV[dim-1][j] << endl;

                // Writing in the recycling file
                for (int k=0; k<n_params; k++){
                    fileRecy << data[ndata+k][j] << ";";
                }
                fileRecy << model[j] << ";";
                for (int k=0; k<diversity-1; k++){
                    fileRecy << recy[k][j] << ";";
                }
                fileRecy << recy[diversity-1][j] << endl;

                // Writing in the recycling CV file
                for (int k=0; k<n_params; k++){
                    fileRecyCV << data[ndata+k][j] << ";";
                }
                fileRecyCV << model[j] << ";";
                for (int k=0; k<diversity-1; k++){
                    fileRecyCV << recyCV[k][j] << ";";
                }
                fileRecyCV << recyCV[diversity-1][j] << endl;

                // Writing in the detritus file
                for (int k=0; k<n_params; k++){
                    fileDetritus << data[ndata+k][j] << ";";
                }
                fileDetritus << model[j] << ";";
                for (int k=0; k<diversity-1; k++){
                    fileDetritus << detritus[k][j] << ";";
                }
                fileDetritus << detritus[diversity-1][j] << endl;

                // Writing in the detritus CV file
                for (int k=0; k<n_params; k++){
                    fileDetritusCV << data[ndata+k][j] << ";";
                }
                fileDetritusCV << model[j] << ";";
                for (int k=0; k<diversity-1; k++){
                    fileDetritusCV << detritusCV[k][j] << ";";
                }
                fileDetritusCV << detritusCV[diversity-1][j] << endl;

                // Writing in the TL file
                for (int k=0; k<n_params; k++){
                    fileTL << data[ndata+k][j] << ";";
                }
                fileTL << model[j] << ";";
                for (int k=0; k<diversity-1; k++){
                    fileTL << TL[k][j] << ";";
                }
                fileTL << TL[diversity-1][j] << endl;

                // Writing in the extinction time file
                for (int k=0; k<n_params; k++){
                    filetExt << data[ndata+k][j] << ";";
                }
                filetExt << model[j] << ";";
                for (int k=0; k<diversity-1; k++){
                    filetExt << tExt[k][j] << ";";
                }
                filetExt << tExt[diversity-1][j] << endl;
            }
        }
        // Delete the result tables
        delete_table(data,ndata+n_params);
        delete_table(biomass,dim);
        delete_table(biomassCV,dim);
        delete_table(recy,diversity);
        delete_table(recyCV,diversity);
        delete_table(detritus,diversity);
        delete_table(detritusCV,diversity);
        delete_table(TL,diversity);
    }

    // Delete the table of parameters
    delete_table(params,n_simu);
    delete_table(pop,1);

    cout << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;

    return 0;
}
