# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")
options(scipen=999)
# import packages
library(lattice)
library(mice)
library(tidyverse)
source('FUNCTIONS.R')

n.imp = 30 # between-variance
nr_boots = 50 #within-variance
# simulate study

study.results = simulation.study(repetitions = 2, 
                                 n.imp = n.imp,
                                 nr_boots = nr_boots, 
                                 mech = 'MAR', 
                                 var = 'covariates')


study.results
capture.output(study.results, file = " ")
# calculate convergence rate
convergence.rate(study.results)

# visualize results
# Opening the graphical device
pdf("plots_MAR_covariates_1000.pdf")
# Creating a plot
plots = forestplot(study.results)
# Closing the graphical device
dev.off() 

## MNAR on outcome
pattern_list = list(c(1,1,0))
prob = .5
mech_list = 'MNAR'


study.results.2 = system.time(simulation.study.loop(repetitions = 5, 
                                                    n.imp = n.imp,
                                                    nr_boots = nr_boots, 
                                                    mech = 'MAR', 
                                                    var = 'outcome'))
