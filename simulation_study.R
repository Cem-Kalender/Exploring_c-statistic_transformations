# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")
options(scipen=999)
# import packages
library(lattice)
library(mice)
library(tidyverse)
source('FUNCTIONS.R')

n.imp = 5 # between-variance

# exectute simulation study
res = simulation.study.Altman(repetitions = 10000, 
                       n.imp = 30, 
                       mech = 'MAR', 
                       var = 'covariates')

# display convergence
convergence.rate(res)
res
# visualize results
forestplot(res)
