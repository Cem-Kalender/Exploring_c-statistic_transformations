# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")
# import packages
source('functions_simulation_study.R')
options(scipen=0)

# Simulation study ----

# set parameters accordingly:
## Scenario: normally distributed MAR missing predictors
#noise = c: (0.25 = 0.7, 0.125 = 0.8, 0.05 = 0.9)
n.imp = 10
rep = 1000
noise_params = c(0.25, 0.125, 0.05)

# Execute simulation study
simulation_study = NULL
for (i in 1:length(noise_params)){
  # set seed to ensure reproducibility
  set.seed(999)
  simulation_study[[i]] = simulation.study(repetitions = rep, 
                                  n.imp = n.imp, 
                                  imp.method = c('norm','norm','logreg'),
                                  noise_param = noise_params[i])
}

# store results of simulation study
saveRDS(simulation_study, file="simulation_study.RData")