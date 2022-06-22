# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")

# import packages
library(lattice)
library(mice)
library(tidyverse)
source('functions_simulation_study.R')

## simulate data set----

## STEP BY STEP ----
n.imp = 10
data = simulate.multivariate(noise=0.05)
plot(data$V1, data$V2, col=ifelse(data$Y==1, 'green', 'red'))
summary(data)
# determine 'true' C-statistic first
sm = sample(seq_len(floor(nrow(data)*0.5)), replace = F)
testdata=data[sm,] # TRUE or FALSE?
data= data[-sm,]
fit = glm(factor(Y)~.,family = 'binomial', data=data)
summary(fit)

pred = predict(fit, newdata=testdata, type= 'response')
true.c = pROC::roc(testdata$Y, pred, quiet = T)$auc[[1]]

true.c
# AMPUTE (Check and confirm with indicator method) FIT 2 MODELS
trainset = amputation(data)
trainset$Y = as.factor(trainset$Y)

testset =  amputation(testdata)
testset$Y = factor(testset$Y)


# IMPUTE ---
train.imp = mice(trainset, method = c('norm','norm','logreg'), n.imp)
test.imp =  mice(testset, method = c('norm','norm','logreg'), n.imp)
# train.imp = mice(trainset, method = c('norm','norm','logreg'), n.imp)
# test.imp = mice(testset, method = c('norm','norm','logreg'), n.imp)

# check convergence
conv = convergence(train.imp)
conv

plot(train.imp)
plot(test.imp)
# fit a model on each imputed data set (SHOULD I USE A SHRINKAGE FACTOR?)
trainset = complete(train.imp, 'long')
testset = complete(test.imp, 'long')
# develop model by combining coefficients with Rubin's rules
model.list = lapply(1:train.imp$m, function(x){
    d = subset(trainset, .imp == x)
    plot(d$V1, d$V2, col=ifelse(d$Y==1, 'green', 'red'))
    d = d %>% select(Y,V2,V1)
    fit = glm(factor(Y)~.,family = 'binomial', data=d) # Shrinkage?
})

final.model = summary(mice::pool(model.list,  rule = "rubin1987"))
# Calculate predicted risk of final model on validation data set ----
testlist <- lapply(1:test.imp$m, function(x){subset(testset, .imp == x)})

# prepare data
validation.data <- lapply(testlist, function(x){
    model.matrix(formula(paste0("~", c('V2+V1'))), data=x)})

predictions = lapply(validation.data, function(x){
    pred = sapply(1:nrow(x), function(i){
        1/(1+exp(-x[i,] %*%  final.model[c(1,2,3),2]))
    })
})

perf.on.validation = lapply(1:n.imp, function(x){
    pred = predictions[[x]]
    test = testlist[[x]]
    auc = pROC::roc(test$Y, pred, quiet = T)$auc[[1]]
    se.auc = auctestr::se_auc(auc, n_p = sum(test$Y==1), n_n = sum(test$Y==0))
    auc_and_se = cbind('C'=auc, 'SE'=se.auc)
    return(auc_and_se)
})

# Results
par(mfrow=c(1,1))
results = data.frame(matrix(unlist(perf.on.validation), nrow=length(perf.on.validation), byrow=TRUE))
colnames(results) <- c('C', 'SE')

# Tranform and pool results
logit.ci = logit.pool(results, n.imp)
regular.ci = regular.pool(results, n.imp)
arcsine.ci = arcsine.pool(results, n.imp)