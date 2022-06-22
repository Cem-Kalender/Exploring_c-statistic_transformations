# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")

# import packages
library(lattice)
library(mice)
library(tidyverse)
source('FUNCTIONS.R')
source('VAL_MI_FUNCTIONS.R')

## simulate dataset----

set.seed(5)
#options(scipen=999)


## SETP BY STEP ----
set.seed(1)
testdata = simulate.multivariate(5000, threshold = 0.5, prob = 0.9)

# determine 'true' C-statistic first
true.c = true.c.stat(testdata)

set.seed(1)
missingdata = amputation(testdata, method = 'MAR', on = 'covariates')

## INDICATOR MATHOD TO CHECK MAR ----
testdata$missV1 = factor(ifelse(is.na(missingdata$V1), 1, 0))
testdata$missV2 = factor(ifelse(is.na(missingdata$V2), 1, 0))

# Fit logistic regression with indicator variable as outcome

glm(factor(missV1)~V1+V2+V3, family = 'binomial', data = testdata) %>% summary()
glm(factor(missV2)~V1+V2+V3, family = 'binomial', data = testdata) %>% summary()

# YES IT WORKS!


# Visualize
ggpubr::ggarrange(
    testdata %>% ggplot(aes(V1,V2, color = factor(V3))) + geom_point() + geom_smooth(method = 'lm'),
    na.omit(missingdata) %>% ggplot(aes(V1,V2, color = factor(V3))) + geom_point() + geom_smooth(method = 'lm')
)

# Impute ---
n.imp=5
missingdata = amputation(testdata, method = 'MAR',on = 'covariates')
missingdata$V3 = as.factor(missingdata$V3)
smp=sample(seq_len(nrow(missingdata)), size = floor(0.5*nrow(missingdata)), replace=T)
train = imputation(missingdata[smp,], n.imp)
test = imputation(missingdata[-smp,], n.imp)
s = sample(seq_len(nrow(train)))


input = list(train, test)

i = processor(input)

sapply(1:n.imp, function(x){
    traindata = i[[1]][[x]][s,]
    model = glm(factor(V3)~V1+V2, family = 'binomial', data = traindata)
    aucs = sapply(1:n.imp, function(z){
        testdata = i[[2]][[z]][s,]
        pred = predict(model, newdata=testdata, type = 'response')
        auc = pROC::roc(testdata$V3, pred, quiet = T)$auc[[1]]
        return(auc)
    })
    return(aucs)
})


input

res = boot::boot(input, statistic = simps, R=100)
mean(res$t)
sd(res$t)

simps = function(missingdata, formula, n.imp){
    # sample with replacement
    smp = sample(seq_len(nrow(missingdata)), replace = T)
    train = missingdata[smp,]

    
    trainlist = impute.missing.values(train, n.imp)
    
    # create model on every list element
    model.list = lapply(trainlist, function(x){
        model = glm(formula, family = 'binomial', data = x)
    })
    
    all.c = sapply(1:n.imp, function(x){
        model = model.list[[x]]
        aucs = sapply(1:n.imp, function(z){
            data = trainlist[[z]]
            pred = predict(model, newdata=data, type = 'response')
            auc = pROC::roc(data$V3, pred, quiet = T)$auc[[1]]
            return(auc)  
        })
        return(aucs)
    })
    c = mean(all.c)
    return(c)
}
simps(missingdata,factor(V3)~V1+V2, 10)

res = boot::boot(missingdata, statistic = simps, formula =factor(V3)~V1+V2, R=100)
results = cbind(C.statistic= mean(res$t), SE=sd(res$t)) %>% data.frame()
results

# forest plot
dist = metamisc::ccalc(cstat = res$t, cstat.se = results$SE)
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


## CHECK ACCURACY OF METHODS ----

# Forest plot
par(mfrow=c(1,1))
metafor::forest(slab = c('Logit', 'Arcsine', 'Regular', 'True'),
                x = combined$C, 
                ci.lb = combined$`95%lb`,
                ci.ub = combined$`95%ub`,
                refline = combined$C[5], xlab = "C-statistic")

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

