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

## SETP BY STEP ----
set.seed(1)
testdata = simulate.multivariate()

# determine 'true' C-statistic first
true.C = boot::boot(test.data,func, R=100)
true.C.CI = boot::boot.ci(true.C, type = 'perc')$percent[c(4,5)]
true = c(true.C.CI[1], true.C$t0, true.C.CI[2])


set.seed(1)
missingdata = amputation(testdata)


# Impute

n.imp = 10
missingdata$V3 = as.factor(missingdata$V3)
set.seed(1)
imp = mice(missingdata, n.imp, method = c('pmm','pmm','logreg'), seed = 5)

# apply bootstrap function and summarize results ----
results = c.stat.bootstrap(imp, nr_boots = 50)
results

# forest plot
dist = metamisc::ccalc(cstat = results$C.statistic, cstat.se = results$SE, N = 20)
plot(dist)


## APPLY TRANSFORMATIONS ----

## LOGIT TRANSFORMATION ----
logit.C = logit.transform.pool(results, n.imp)
arcsin.C = arcsine.pool(results, n.imp)
regular = regular.pool(results, n.imp)

## Compare all transformed methods
combined = rbind(logit.C, arcsin.C, regular, true) %>% data.frame()
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
metafor::forest(slab = c('Logit', 'Arcsine', 'Regular', 'True'),
                x = combined$C, 
                ci.lb = combined$`95%lb`,
                ci.ub = combined$`95%ub`,
                refline = combined$C[5], xlab = "C-statistic")

