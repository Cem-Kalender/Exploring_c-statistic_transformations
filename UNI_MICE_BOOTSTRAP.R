# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")

# import packages
library(lattice)
library(mice)
library(tidyverse)

## simulate dataset----

set.seed(5)
#options(scipen=999)

## read data
test.data = read.csv('complete_data.csv', header = T, sep = ',')
missingdata = read.csv('univariate_missingdata.csv', header = T, sep = ',')
is.na(missingdata) %>% sum()


## FUNCTIONS ----

#define function for bootstrapping
func = function(data,ind){
    fit = glm(factor(V3) ~ V2 + V1, data = data[ind,], family = "binomial")
    preds = predict(fit, type = 'response', newdata = test.data)
    pROC::auc(pROC::roc(test.data$V3, preds))
    #DescTools::Cstat(fit)
}


## Function for SE pooling
pooled_se <- function(est, se, nimp){
    w_auc <- 
        mean(se^2) # within variance
    b_auc <-
        var(est) # between variance
    tv_auc <-
        w_auc + (1 + (1/nimp)) * b_auc # total variance
    se_total <-
        sqrt(tv_auc)
    r <- (1 + 1 / nimp) * (b_auc / w_auc)
    v <- (nimp - 1) * (1 + (1/r))^2
    t <- qt(0.975, v)
    res <- c(se_total, t)
    return(res)
}



## MI ----

# determine 'true' C-statistic first
true.model = glm(factor(V3) ~ V2 + V1, data = testdata, family = "binomial")
true.model %>% summary()

true.C = boot::boot(test.data,func, R=100)
true.C.CI = boot::boot.ci(true.C, type = 'perc')$percent[c(4,5)]
true = c(true.C$t0, true.C.CI[1], true.C.CI[2])
true

# fit a model based on incomplete data and inspect
bias.model = glm(factor(V3) ~ V2 + V1, data = missingdata, family = "binomial")
bias.model %>% summary()

biased.C = boot::boot(missingdata, func, R=100)
biased.C.CI = boot::boot.ci(biased.C, type = 'perc')$percent[c(4,5)]
biased = c(biased.C$t0, biased.C.CI[1], biased.C.CI[2])
biased

# Create empty matrix to store bootstrapped C-statistic estimates
results = matrix(nrow = n.imp, ncol = 2) %>% data.frame()
colnames(results) = c('C.statistic', 'SE')

# Impute
set.seed(5)
missingdata$V3 = as.factor(missingdata$V3)
n.imp = 20
imp = mice(missingdata, n.imp, method = c('norm.nob','norm.nob', 'logreg'), seed = 5) # stochastic
imp

# stack
imp_tot <- complete(imp, "stacked")
start = seq(1, 160000, by=8000)
end = seq(8000, 160000, by=8000)
start
end

# Multiple Imputation
for (i in 1:n.imp){
    
    data.copy = imp_tot[start[i]:end[i],]
    
    ## NON-PARAMETRIC BOOTSTRAP for determining within-variance of C-statistic
    
    ## create bootstrap
    res = boot::boot(data.copy, func, R=100)
    # store sample statistic and corresponding error (SD=SE)
    results[i, 1] = res$t0
    results[i, 2] = sd(res$t)
    
    # Use imputed data to fit a prediction model (and estimate performance on test data)
    #pred.model = glm(factor(V3, c(0,1))~V2+V1, family = 'binomial', data.copy)
    
    # predict
    #preds = predict(pred.model, type = 'response', newdata = test.data)
    
    # calculate C-statistic and insert
    #ROC = pROC::roc(test.data$V3, preds)
    #results[i, 3] <- pROC::auc(ROC)
    
}

# inspect results----
summary(results)
results

## VISUALIZE TRANSFORMATIONS ----

par(mfrow=c(1,3))
# plot histogram
hist(results$C.statistic, freq = F, breaks = 10, main = 'Distribution of C-statistic')
# plot densityplot
dens <- density(results$C.statistic)
lines(dens, col = 'blue')

#logit transformation of the c-statistics
logit =  log(results/(1 - results))
# plot histogram
hist(logit$C.statistic, freq = F, breaks = 10, main = 'Distribution of logit transformed C-statistic')
# plot densityplot
dens <- density(logit$C.statistic)
lines(dens, col = 'red')

# ARCSINE
arcsin =  asin(sqrt(results$C.statistic))
# plot histogram
hist(arcsin, freq = F, breaks = 10, main = 'Distribution of arcsine transformed C-statistic')
# plot densityplot
arcsin.dens <- density(arcsin)
lines(arcsin.dens, col = 'green')


## LOGIT TRANSFORMATION ----

logit.pooled.se = pooled_se(results$C.statistic, results$SE, n.imp)

logit.ub = exp(mean(logit$C.statistic) + (logit.pooled.se[2]*logit.pooled.se[1])) /
    (1 + exp(mean(logit$C.statistic) + (logit.pooled.se[2]*logit.pooled.se[1])))

logit.lb = exp(mean(logit$C.statistic) - (logit.pooled.se[2]*logit.pooled.se[1])) /
    (1 + exp(mean(logit$C.statistic) - (logit.pooled.se[2]*logit.pooled.se[1])))

# Inspect

logit.C = c(Pooled.C.statistic = exp(mean(logit$C.statistic))/(1 + exp(mean(logit$C.statistic))), 
            '95%.Low' =logit.lb,
            '95%.Up' =logit.ub)


## ARCSINE TRANSFORMATION ----

# Pooled SE
arcsin.pooled.se = pooled_se(results$C.statistic, results$SE, n.imp)
# Upper bound
arcsin.ub = sin(mean(arcsin))**2 + (arcsin.pooled.se[2]*arcsin.pooled.se[1])/sqrt(n.imp)
# lower bound
arcsin.lb = sin(mean(arcsin))**2 - (arcsin.pooled.se[2]*arcsin.pooled.se[1])/sqrt(n.imp)

# Inspect
arcsin.C = c(Pooled.C.statistic = sin(mean(arcsin))**2, 
             '95%.Low' = arcsin.lb,
             '95%.Up' = arcsin.ub)


## REGULAR POOLING ----

## pooled C-statistic estimates is just the average:
pooled.se = pooled_se(results$C.statistic, results$SE, n.imp)
# Inspect
C = c(Pooled.C.statistic = mean(results$C.statistic), 
      '95%.Low' = mean(results$C.statistic) - (1.96*pooled.se[1]/sqrt(n.imp)),
      '95%.Up' = mean(results$C.statistic) + (1.96*pooled.se[1]/sqrt(n.imp)))

## Compare all transformed methods
combined = rbind(logit.C, arcsin.C, C, true) %>% data.frame()
combined

## CHECK ACCURACY OF METHODS ----

# Forest plot
par(mfrow=c(1,1))
metafor::forest(slab = c('Logit', 'Arcsine', 'regular', 'True'),
                x = combined$Pooled.C.statistic, 
                ci.lb = combined$X95..Low,
                ci.ub = combined$X95..Up,
                refline = combined$Pooled.C.statistic[4], xlab = "C-statistic")

