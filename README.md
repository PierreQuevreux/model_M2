# model_M2
Models used in the article "Interplay between the paradox of enrichment and nutrient cycling in food webs"

Program coded in C++ and using the GSL ODE solver (v 1.13).

# Foodchain
Simulation of a food chain. The program records the dynamics of the biomasses and nutrients.

# Foodweb
Simulation of a complex food web with up to 50 interacting species. The global flows of matter and the average biomasses are recorded. Time dynamics of species are also recorded. The code uses the openMP library for parallel calculation.

# Figures.R
R code to produce the figures of the article

Each program stores the data in a "data" or "results" folder containing a README.txt detailing the recorded variables and parameters.
