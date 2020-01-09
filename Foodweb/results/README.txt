Files containing the data

### Files ###

parameters_generator.R : create the parameters_N.txt and parameters_data_N.txt files with parameters values
aggregation_data.R : merge the output files of several simulations

###########
data_N.txt : File containing the general data of the simulation

NbSpeciesInit : initial number of species
NbPPInit : initial number of primary producers
NbSpeciesFinal : final number of surviving species
NbPPFinal : final number of surviving primary producers
TLmax : maximum trophic level in the food web
connectance : connectance
totalBiomass : average total biomass of all species
PPbiomass : average primary producers total biomass
SPbiomass : average consumer total biomass
totalProd : average total biomass production
PPprod : average primary production
SPprod : average secondary production
Irecy : average quantity of recycled nutrients
RecyPP : average quantity of nutrients recycled by primary producers
RecySP : average quantity of nutrients recycled by consumers
RecyDir : average quantity of directly recycled nutrients
RecyInd : average quantity of indirectly recycled nutrients
totalMort : average total biomass loss due to density dependent mortality
PPmort : average total biomass loss due to density dependent mortality of primary producers
SPmort : average total biomass loss due to density dependent mortality of consumers
totalBiomassCV : coefficient of variation of total biomass
PPbiomassCV : coefficient of variation of primary producer total biomass
SPbiomassCV : coefficient of variation of consumer total biomass
totalProdCV : coefficient of variation of total biomass production
PPprodCV : coefficient of variation of primary prodcution
SPprodCV : coefficient of variation of secondary production
IrecyCV : coefficient of variation of quantity of recycled nutrients
RecyPPcv : coefficient of variation of primary production
RecySPcv : coefficient of variation of secondary production
RecyDirCV : coefficient of variation of the quantity of directly recycled nutrients
RecyIndCV : coefficient of variation of the quantity of indirectly recycled nutrients
totalMortCV : coefficient of variation of the total biomass loss due to density dependent mortality
PPmortCV : coefficient of variation of the total biomass loss due to density dependent mortality of primary producers
SPmortCV : coefficient of variation of the total biomass loss due to density dependent mortality of consumers
nbSimu : ID of the simulation
seed : seed of the random number generator
I : nutrient input
delta : fraction of direct recycling
d : decomposition rate
model : NC, C or SC

###########
bodymass_N.txt : File containing the bady mass of each species

nbSimu : ID of the simulation
xi : species i

###########
biomass_N.txt : File containing the average biomass of each species

nbSimu : ID of the simulation
seed : seed of the random number generator
I : nutrient input
delta : fraction of direct recycling
d : decomposition rate
model : NC, C or SC
N : mineral nutrients
D : detritus
xi : species i

###########
biomassCV_N.txt : File containing the biomass CV of each species

nbSimu : ID of the simulation
seed : seed of the random number generator
I : nutrient input
delta : fraction of direct recycling
d : decomposition rate
model : NC, C or SC
N : mineral nutrients
D : detritus
xi : species i

###########
recy_N.txt : File containing the average quantity of nutrient directly recycled by each species

nbSimu : ID of the simulation
seed : seed of the random number generator
I : nutrient input
delta : fraction of direct recycling
d : decomposition rate
model : NC, C or SC
xi : species i

###########
recyCV_N.txt : File containing the CV of the quantity of nutrient directly recycled by each species

nbSimu : ID of the simulation
seed : seed of the random number generator
I : nutrient input
delta : fraction of direct recycling
d : decomposition rate
model : NC, C or SC
xi : species i

###########
detritus_N.txt : File containing the average quantity of detritus released by each species

nbSimu : ID of the simulation
seed : seed of the random number generator
I : nutrient input
delta : fraction of direct recycling
d : decomposition rate
model : NC, C or SC
xi : species i

###########
detritusCV_N.txt : File containing the CV of the quantity of detritus released by each species

nbSimu : ID of the simulation
seed : seed of the random number generator
I : nutrient input
delta : fraction of direct recycling
d : decomposition rate
model : NC, C or SC
xi : species i

###########
TL_N.txt : File containing the trophic level of each species

nbSimu : ID of the simulation
seed : seed of the random number generator
I : nutrient input
delta : fraction of direct recycling
d : decomposition rate
model : NC, C or SC
xi : species i

###########
tExtinction_N.txt : File containing the extinction time of each species

nbSimu : ID of the simulation
seed : seed of the random number generator
I : nutrient input
delta : fraction of direct recycling
d : decomposition rate
model : NC, C or SC
xi : species i

###########
time_series_data_N.txt : File containing the time series of the main variables

nbSimu : ID of the simulation
seed : seed of the random number generator
I : nutrient input
delta : fraction of direct recycling
d : decomposition rate
model : NC, C or SC
t : time (years)
PP : primary production
SP : secondary production
Irecy : total quantity of recycled nutrients
RecyPP : total quantity of nutrients recycled by primary producres
RecySP : total quantity of nutrients recycled by consumers
RecyInd : total indirect recycling
RecyDir : total direct recycling
totalBiomass : total biomass
PPbiomass : primary producer biomass
SPbiomass : consumer biomass

###########
time_series_biomass_N.txt : File containing the time series of the biomass of each species

nbSimu : ID of the simulation
seed : seed of the random number generator
I : nutrient input
delta : fraction of direct recycling
d : decomposition rate
model : NC, C or SC
t : time (years)
N : mineral nutrients
D : detritus
Xi : species i

###########
time_series_recy_N.txt : File containing the time series of the nutrient directly recycled by each species

nbSimu : ID of the simulation
seed : seed of the random number generator
I : nutrient input
delta : fraction of direct recycling
d : decomposition rate
model : NC, C or SC
t : time (years)
Xi : species i

###########
time_series_detritus_N.txt : File containing the time series of the quantity of detritus released by each species

nbSimu : ID of the simulation
seed : seed of the random number generator
I : nutrient input
delta : fraction of direct recycling
d : decomposition rate
model : NC, C or SC
t : time (years)
Xi : species i
