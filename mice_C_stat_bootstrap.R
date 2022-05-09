# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")

# import packages
library(lattice)
library(mice)
library(tidyverse)
source('FUNCTIONS.R')

## simulate dataset----

set.seed(5)
#options(scipen=999)

## read data
test.data = read.csv('complete_data.csv', header = T, sep = ',')
missingdata = read.csv('missingdata.csv', header = T, sep = ',')
# inspect data
summary(missingdata)

# Visualize
ggpubr::ggarrange(
    test.data %>% ggplot(aes(V1,V2, color = factor(V3))) + geom_point() + geom_smooth(method = 'lm'),
    na.omit(missingdata) %>% ggplot(aes(V1,V2, color = factor(V3))) + geom_point() + geom_smooth(method = 'lm')
)

# determine 'true' C-statistic first
true.model = glm(factor(V3) ~ V2 + V1, data = testdata, family = "binomial")
true.model %>% summary()

true.C = boot::boot(test.data,func, R=100)
true.C.CI = boot::boot.ci(true.C, type = 'perc')$percent[c(4,5)]
true = c(true.C.CI[1], true.C$t0, true.C.CI[2])
true

# fit a model based on incomplete data and inspect
bias.model = glm(factor(V3) ~ V2 + V1, data = missingdata, family = "binomial")
bias.model %>% summary()

biased.C = boot::boot(missingdata, func, R=100)
biased.C.CI = boot::boot.ci(biased.C, type = 'perc')$percent[c(4,5)]
biased = c(biased.C.CI[1],biased.C$t0, biased.C.CI[2])
biased

# Impute
set.seed(5)
n.imp = 20
missingdata$V3 = as.factor(missingdata$V3)
imp = mice(missingdata, n.imp, method = c('norm.nob','norm.nob', 'logreg'), seed = 5) # stochastic

# fmri----
imp_tot$Impnr = factor(rep(1:20, each = 10000))

perf <- psfmi::pool_performance(data=imp_tot, nimp=n.imp, impvar="Impnr", 
                         formula = V3 ~ V1 + V2, 
                         cal.plot=TRUE, plot.method="mean", 
                         groups_cal=20, model_type="binomial")

perf$ROC_pooled

# apply bootstrap function and summarize results
results = C_stat.bootstrap(imp, nr_boots = 100)
summary(results)

## APPLY TRANSFORMATIONS ----

## MEDIAN 
med = c(median(results$C.statistic), IQR(results$C.statistic))
median = c('95%.Low' = med[1] - 1.96*med[2]/sqrt(n.imp),Pool3e.C.statistic = med[1], '95%.Up' = med[1] + 1.96*med[2]/sqrt(n.imp))


## LOGIT TRANSFORMATION ----
logit =  log(results/(1 - results))
est_c_se_log = results$SE / (results$C.statistic * (1-results$C.statistic))

logit.pooled.se=pooled_se(logit$C.statistic, est_c_se_log, n.imp)

logit.ub = exp(mean(logit$C.statistic) + (logit.pooled.se[2]*logit.pooled.se[1])) /
    (1 + exp(mean(logit$C.statistic) + (logit.pooled.se[2]*logit.pooled.se[1])))

logit.lb = exp(mean(logit$C.statistic) - (logit.pooled.se[2]*logit.pooled.se[1])) /
    (1 + exp(mean(logit$C.statistic) - (logit.pooled.se[2]*logit.pooled.se[1])))

# Inspect

logit.C = c('95%.Low' = logit.lb,
            Pooled.C.statistic = exp(mean(logit$C.statistic))/(1 + exp(mean(logit$C.statistic))),
            '95%.Up' =logit.ub)
logit.C

## ARCSINE TRANSFORMATION 

# Pooled SE
arcsin =  asin(sqrt(results$C.statistic))
est_c_se_arc = results$SE/(arcsin * (1-arcsin))

arcsin.pooled.se = pooled_se(arcsin, est_c_se_arc, n.imp)
# Upper bound
arcsin.ub = sin(mean(arcsin))**2 + (arcsin.pooled.se[2]*arcsin.pooled.se[1])/sqrt(n.imp)
# lower bound
arcsin.lb = sin(mean(arcsin))**2 - (arcsin.pooled.se[2]*arcsin.pooled.se[1])/sqrt(n.imp)

# Inspect
arcsin.C = c('95%.Low' = arcsin.lb,
             Pooled.C.statistic = sin(mean(arcsin))**2,
             '95%.Up' = arcsin.ub)

arcsin.C

## REGULAR POOLING 

## pooled C-statistic estimates is just the average:
pooled.se = pooled_se(results$C.statistic, 
                      results$SE/(results$C.statistic*(1-results$C.statistic)), 
                      n.imp)
# Inspect
regular = c('95%.Low' = mean(results$C.statistic) - (pooled.se[2]*pooled.se[1]/sqrt(n.imp)),
            Pooled.C.statistic = mean(results$C.statistic),
            '95%.Up' = mean(results$C.statistic) + (pooled.se[2]*pooled.se[1]/sqrt(n.imp)))

## Compare all transformed methods
combined = rbind(logit.C, arcsin.C, regular, median, true) %>% data.frame()
colnames(combined) = c('95%lb', 'C', '95%ub')
combined

## VISUALIZE TRANSFORMATIONS ----

par(mfrow=c(1,3))

# plot histogram
hist(results$C.statistic, freq = F, breaks = 10, main = 'Distribution of C-statistic')
# plot densityplot
dens <- density(results$C.statistic)
lines(dens, col = 'blue')

# plot histogram
hist(logit$C.statistic, freq = F, breaks = 10, main = 'Distribution of logit transformed C-statistic')
# plot densityplot
dens <- density(logit$C.statistic)
lines(dens, col = 'red')

# plot histogram
hist(arcsin, freq = F, breaks = 10, main = 'Distribution of arcsine transformed C-statistic')
# plot densityplot
arcsin.dens <- density(arcsin)
lines(arcsin.dens, col = 'green')


## CHECK ACCURACY OF METHODS ----

# Forest plot
par(mfrow=c(1,1))
metafor::forest(slab = c('Logit', 'Arcsine', 'Regular', 'Median', 'True'),
                x = combined$C, 
                ci.lb = combined$`95%lb`,
                ci.ub = combined$`95%ub`,
                refline = combined$C[5], xlab = "C-statistic")

