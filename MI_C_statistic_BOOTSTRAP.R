# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")

# import packages
library(lattice)
library(mice)
library(tidyverse)

## simulate dataset----

set.seed(5)
#options(scipen=999)

## read amputed data
missingdata = read.csv('missingdata.csv', header = T, sep = ',')

## MI ----

N = sum(is.na(missingdata$V1))
miss = missingdata[is.na(missingdata$V1),]
n.imp = 20
results = matrix(nrow = n.imp, ncol = 3) %>% data.frame()
colnames(results) = c('C.statistic', '95%CI.lower', '95%CI.upper')

#define function for bootstrapping
func = function(data,ind){
    fit = glm(factor(V3) ~ V2 + V1, data = data[ind,], family = "binomial")
    DescTools::Cstat(fit)
}

 # Multiple Imputation
for (i in 1:n.imp){
    # copy missing data
    data.copy = missingdata
    
    # create imputation model
    imp.model = lm(V1~factor(V3, c(0,1))+V2, data=data.copy)
    
    #drawing coefficients from a multivariate Student-t distribution
    betas <- mnormt::rmt(N, mean=imp.model$coef, S=vcov(imp.model), df=imp.model$df.residual)
    
    # retrieve linear predictors
    predictors = model.matrix(formula(paste("~", 'factor(V3, c(0,1))+V2')), data=miss)
    
    # empty imputations object to store imputed values
    imputations = NULL
    
    # matrix multiply predictors with betas for j subjects
    # + add noise by drawing from residual distribution
    for (j in 1:N){
        imp.value = (predictors[j,] %*% betas[j,]) + sample(imp.model$residuals, 1)
        #imp.prob.value = 1/(1+exp(-predictors[j,] %*% betas[j,]))
        imputations = append(imputations, imp.value)
        #probs = append(probs, imp.prob.value)
    }
    
    #imputations = rbinom(n = N, size = 1, prob = probs)
    # Impute data
    data.copy$V1[is.na(data.copy$V1)] = imputations
    
    ## NON-PARAMETRIC BOOTSTRAP for C-statistic
    
    ## create bootstrap
    res = boot::boot(data.copy, func, R=100)
    # Percentiles
    C = boot::boot.ci(res,type="perc")
   
    # take percentile interval of C-statistic
    results[i, c(2,3)] = C$percent[c(4,5)]
    results[i, 1] = C$t0
    
}

# inspect results----
summary(results)

## VISUALIZE TRANSFORMATIONS ----

par(mfrow=c(1,3))
# plot histogram
hist(results$C.statistic, freq = F, breaks = 10, main = 'Distribution of C-statistic')
# plot densityplot
dens <- density(results$C.statistic)
lines(dens, col = 'blue')

#logit transformation of the c-statistics
logit.C =  log(results/(1 - results))
# plot histogram
hist(logit.C$C.statistic, freq = F, breaks = 10, main = 'Distribution of logit transformed C-statistic')
# plot densityplot
dens <- density(logit.C$C.statistic)
lines(dens, col = 'red')

#logit transformation of the c-statistics
arcsine.C =  asin(sqrt(results$C.statistic))

# plot histogram
hist(arcsine.C, freq = F, breaks = 10, main = 'Distribution of arcsine transformed C-statistic')
# plot densityplot
arcsine.dens <- density(arcsine.C)
lines(arcsine.dens, col = 'green')


## apply Rubin's rules on logit transformed C-statistic ----

## pooled C-statistic estimates is just the average:
beta = logit.C$C.statistic %>% mean()

# retrieve standard error of bootstrapped CI estimate 
SE = sqrt(nrow(missingdata + 1)) * (logit.C[,3] - logit.C[,1])/1.96
SE

# The SE for the pooled C-statistic estimate is retrieved using Rubin's rules
Qbar <- SE %>% mean()
U <- sum(SE**2)/n.imp
B <- sum((logit.C$C.statistic - Qbar)**2)/(n.imp-1)
se.beta <- sqrt(U + (1+1/n.imp)*B)
# Inspect
C = c(Pooled.C.statistic =exp(beta)/(1+exp(beta)), SE =exp(se.beta)/(1+exp(se.beta)))
C
