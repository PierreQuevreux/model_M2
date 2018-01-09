# model_M2
Models used in the article "Impact of nutrient cycling on food web stability"

Program coded in C++ and using the GSL ODE solver (v 1.13).

# Foodchain
Simulation of a tri-trophic food chain with a primary producer, a herbivore and a carnivore. The program records the dynamics of the biomasses and nutrients. The data for a diagram of bifurcation are also produced (calculations of local extrema along parameter gradients).

# Foodweb
Simulation of a complex food web with up to 50 interactin species. The global flows of matter and the average biomasses are recorded. Time dynamics of species are also recorded.

# Foodweb_simulation
Simulations testing the effect of nutrient cycling along an enrichment gradient. Only the general data such as average species biomasses, coefficients of variation or global flows are recorded.

Each program stores the data in a "data" or "results" folder containing a README.txt detailing the recorded variables and parameters.
