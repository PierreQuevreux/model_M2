Files containing the data

###########
data.txt : File containing the dynamics of biomass
t : time (year)
N : mineral nutrients
D : detritus
PP : primary production
SP : secondary production
PR : primary recycling
SR : secondary recycling
TR : total recycling
Xi : biomass of species i

###########
recy.txt : File containing the dynamics of recycling
t : time (year)
N : mineral nutrients
D : detritus
PP : primary production
SP : secondary production
PR : primary recycling
SR : secondary recycling
TR : total recycling
Xi : quantity of nutrients recycled by species i

###########
pref.txt : File containing the dynamics of preferences of species for their prey
# first line
t : time (year)
i : number of the consumer (x number of prey)
# first line
t : time (year)
j : number of the prey

###########
general.txt : File containing the generalinformation on species
Species : ID of the species
Type : primary producer or consumer
Mass : body mass
TL : trophic level
Average biomass : average biomass of species
Extinction time : >0 if the species got extinct
Cv : coefficient of variation of species biomass
