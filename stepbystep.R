# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")

# import packages
library(lattice)
library(mice)
library(tidyverse)
source('FUNCTIONS.R')

## simulate data set----

## STEP BY STEP ----
n.imp = 10

data = simulate.multivariate(noise = 0.25)
plot(data$V1, data$V2, col=ifelse(data$Y==1, 'green', 'red'))
summary(data)
# determine 'true' C-statistic first
sm = sample(seq_len(floor(nrow(data)*0.5)), replace = F)
testdata=data[sm,] # TRUE or FALSE?
data= data[-sm,]
fit = glm(factor(Y)~.,family = 'binomial', data=data)


pred = predict(fit, newdata=testdata, type= 'response')
true.c = pROC::roc(testdata$Y, pred, quiet = T)$auc[[1]]

# AMPUTE (Check and confirm with indicator method) FIT 2 MODELS
trainset = amputation(data, method = 'MAR', on = 'X&Y')
trainset$Y = as.factor(trainset$Y)

testset =  amputation(testdata, method = 'MAR', on = 'X&Y')
testset$Y = factor(testset$Y)
# Split
# smp = sample(seq_len(floor(nrow(missingdata)*0.5)), replace = F) # TRUE or FALSE?
# trainset = missingdata[smp,]
# testset =  missingdata[-smp,]
# IMPUTE ---
train.imp = imputation(trainset, imp_method = c('norm','norm','logreg'), n.imp)
test.imp = imputation(testset, imp_method = c('norm','norm','logreg'), n.imp)

plot(train.imp)
plot(test.imp)
# fit a model on each imputed data set (SHOULD I USE A SHRINKAGE FACTOR?)
trainset = complete(train.imp, 'long')
testset = complete(test.imp, 'long')
# develop model by combining coefficients with Rubin's rules
model.list = lapply(1:train.imp$m, function(x){
    d = subset(trainset, .imp == x)
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
results = data.frame(matrix(unlist(perf.on.validation), nrow=length(perf.on.validation), byrow=TRUE))
colnames(results) <- c('C', 'SE')
results

# Tranform and pool results
logit.ci = logit.transform.pool(results, n.imp)
regular.ci = regular.pool(results, n.imp)
arcsine.ci = arcsine.pool(results, n.imp)
true.ci = rep(true.c, 3)

# Combine
study.result = rbind(logit.ci, regular.ci, arcsine.ci, true.ci)
study.result

