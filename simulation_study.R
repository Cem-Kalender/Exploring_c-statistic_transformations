# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")
# import packages
source('functions_simulation_study.R')

## Scenario: normally distributed MAR missing predictors
#noise = c: (0.425 = 0.6, 0.25 = 0.7, 0.125 = 0.8, 0.05 = 0.9)
# m appears to affect convergence rate positively
n.imp = 10
rep = 1000
noise_params = c(0.425, 0.25, 0.125, 0.05)

# execute simulation study
simulation_study = NULL
for (i in 1:length(noise_params)){
    rex = simulation.study(repetitions = rep, 
                                  n.imp = n.imp, 
                                  mech = 'MAR',
                                  imp.method = c('norm','norm','logreg'),
                                  var = 'X&Y',
                                  noise_param = noise_params[i])
    simulation_study[[i]]=list('coverage'=coverage.rate(rex), 'bias'=bias(rex))
    }

# What is the coverage rate
coverage_df = sapply(1:4, function(i){simulation_study[[i]]$coverage}) %>% data.frame()
colnames(coverage_df) = c('true.c.0.6', 'true.c.0.7', 'true.c.0.8', 'true.c.0.9')
coverage_df

# What is the bias?
bias_df = sapply(1:4, function(i){simulation_study[[i]]$bias}) %>% data.frame()
colnames(bias_df) = c('true.c.0.6', 'true.c.0.7', 'true.c.0.8', 'true.c.0.9')
apply(bias_df, 2, mean)