Files containing the data

### Files ###

parameters_generator.R : create the parameters_N.txt and parameters_data_N.txt files with parameters values
aggregation_data.R : merge the output files of several simulations
script_simu : bash script to launch the simulations

### Recorded variables ###

simu_ID : ID of the simulation
I : nutrient input
delta : fraction of direct recycling
d : decomposition rate
model : NC, C or SC
t : time (year)
N : mineral nutrients
D : detritus
Xi : biomass of species i
GD_Xi : growth:death ratio of species i
Prod_Xi : biomass production of species i
Recy_Xi : quantity of nutrient directly recycled by species i
IndRecy : nutrient indirectly recycled (detritus decomposition)
SProd : secondary production
SRecy : consumer direct recycling
TProd : total biomass production
TDirRecy : total direct recycling
TRecy : total recycling

### File names ###

time_series_N.txt : time series of the recorded variables
mean_N.txt : mean values of the recorded variables
CV_N.txt : coefficient of variation of the recorded variables
