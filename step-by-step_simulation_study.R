# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")

# import packages
library(lattice)
library(mice)
library(tidyverse)
source('FUNCTIONS.R')
#source('VAL_MI_FUNCTIONS.R')

## simulate data set----

## STEP BY STEP ----
set.seed(1)
testdata = simulate.multivariate(10000, threshold = 0.6, prob = 0.9)
stargazer::stargazer(testdata, type= "html", title= "Summary Statistics", out= "dat.html",
                     covariate.labels = c('A', 'B', 'Y'))
plot(testdata$V1, testdata$V2, col = ifelse(testdata$V3==1,'green', 'red'))

# AMPUTE (Check and confirm with indicator method) FIT 2 MODELS
missingdata = amputation(testdata, method = 'MAR', on = 'covariates')
missingdata$V3 = as.factor(missingdata$V3)

# Visualize 
ggpubr::ggarrange(
    testdata %>% ggplot(aes(V1,V2, color = factor(V3))) + geom_point() + geom_smooth(method = 'lm'),
    na.omit(missingdata) %>% ggplot(aes(V1,V2, color = factor(V3))) + geom_point() + geom_smooth(method = 'lm')
)


# determine 'true' C-statistic first
true.c = true.c.stat(testdata)

# IMPUTE ---

n.imp = 10
smp = sample(seq_len(floor(nrow(missingdata)*0.5)), replace = F) # TRUE or FALSE?
trainset = missingdata[smp,]
testset =  missingdata[-smp,]

train.imp = imputation(trainset, n.imp)
test.imp = imputation(testset, n.imp)

# check for non-convergence
plot(train.imp) # How to automate check for convergence with convergence()?
#convergence(train.imp) How does this function work?
plot(test.imp)

# fit a model on each imputed data set (SHOULD I USE A SHRINKAGE FACTOR?)
trainset = complete(train.imp, 'long')
testset = complete(test.imp, 'long')

# develop model by combining coefficients with Rubin's rules
model.list = lapply(1:train.imp$m, function(x){
    d = subset(trainset, .imp == x)
    d = d %>% select(V3,V2,V1)
    fit = glm(factor(V3)~.,family = 'binomial', data=d) # Shrinkage?
    })

final.model = summary(pool(model.list,  rule = "rubin1987"))
final.model # whats wrong with the p value..?

# CHECK PERFORMANCE ON DEVELOPMENT AND VALIDATION SETS ----
# Calculate predicted risk per model on complete data set
test.on.developmentset = function(model.list){
    c = lapply(model.list, function(x){
        auc.ci = pROC::ci.auc(x$y, x$fitted.values, quiet = T)
        return(auc.ci)
    })
    # reformat
    c = data.frame(matrix(unlist(c), nrow=length(c), byrow=TRUE))
    colnames(c) <- c('95%lb', 'C', '95%ub')
    
    # Tranform and pool results
    logit.ci = logit.transform.pool(c, n.imp)
    regular.ci = regular.pool(c, n.imp)
    
    # Combine
    study.result = rbind(logit.ci, regular.ci)
    return(study.result)
}
test.performance = function(imp, complete.set, final.model, n.imp){
    alist <- lapply(1:imp$m, function(x){subset(complete.set, .imp == x)})
    
    # prepare data
    validation.data <- lapply(alist, function(x){
        model.matrix(formula(paste0("~", c('V2+V1'))), data=x)})
    
    
    predictions = lapply(validation.data, function(x){
        predictions = NULL
        for (i in 1:nrow(x)){
            pred = 1/(1+exp(x[i,] %*%  final.model[c(1,2,3),2]))
            predictions = append(predictions, pred)
        }
        return (predictions)
    })
    
    perf.on.validation = lapply(1:n.imp, function(x){
        pred = predictions[[x]]
        test = alist[[x]]
        auc = pROC::roc(test$V3, pred, quiet = T)$auc[[1]]
        # is the SE unbiased?
        #se.auc = auctestr::se_auc(auc, n_p = sum(test$V3==1), n_n = sum(test$V3==0))
        auc.ci = pROC::ci.auc(test$V3, pred, quiet = T)
        #auc_and_se = cbind('C'=auc, auc.ci)
        return(auc.ci)
    })
    
    # Results
    results = data.frame(matrix(unlist(perf.on.validation), nrow=length(perf.on.validation), byrow=TRUE))
    colnames(results) <- c('95%lb', 'C', '95%ub')
    
    # Tranform and pool results
    logit.ci = logit.transform.pool(results, n.imp)
    regular.ci = regular.pool(results, n.imp)
    
    # Combine
    study.result = rbind(logit.ci, regular.ci)
    
    return(study.result)
    
}

val.on.development = 
    test.on.developmentset(model.list)

val.on.validation = 
    test.performance(test.imp, testset, final.model, n.imp)

# Forest plot
par(mfrow=c(1,1))
metafor::forest(slab = c('Logit transformed', 'Regular pooled'),
                x = val.on.development[,2], 
                ci.lb = val.on.development[,1],
                ci.ub = val.on.development[,3],
                refline = true.c, 
                xlab = "C-statistic",
                main = 'Performance on development set')

metafor::forest(slab = c('Logit transformed', 'Regular pooled'),
                x = val.on.validation[,2], 
                ci.lb = val.on.validation[,1],
                ci.ub = val.on.validation[,3],
                refline = true.c, 
                xlab = "C-statistic",
                main = 'Performance on validation set')


# Calculate predicted risk of final model on validation data set ----
testlist <- lapply(1:test.imp$m, function(x){subset(testset, .imp == x)})
testlist

# prepare data
validation.data <- lapply(testlist, function(x){
    model.matrix(formula(paste0("~", c('V2+V1'))), data=x)})


predictions = lapply(validation.data, function(x){
    predictions = NULL
    for (i in 1:nrow(x)){
        pred = 1/(1+exp(x[i,] %*%  final.model[c(1,2,3),2]))
        predictions = append(predictions, pred)
    }
    return (predictions)
    })

perf.on.validation = lapply(1:n.imp, function(x){
    pred = predictions[[x]]
    test = testlist[[x]]
    auc = pROC::roc(test$V3, pred, quiet = T)$auc[[1]]
    # is the SE unbiased?
    #se.auc = auctestr::se_auc(auc, n_p = sum(test$V3==1), n_n = sum(test$V3==0))
    auc.ci = pROC::ci.auc(test$V3, pred, quiet = T)
    #auc_and_se = cbind('C'=auc, auc.ci)
    return(auc.ci)
})

# Results
results = data.frame(matrix(unlist(perf.on.validation), nrow=length(perf.on.validation), byrow=TRUE))
colnames(results) <- c('95%lb', 'C', '95%ub')
results <- rbind(results, rep(true.c, 3))
# Check distribution
hist(results$C, breaks = 10)
#combined = transform_and_pool(results, true.c, n.imp)

par(mfrow=c(1,1))
metafor::forest(x = results$C, 
                ci.lb = results$`95%lb`,
                ci.ub = results$`95%ub`,
                refline = true.c, 
                xlab = "C-statistic")

# Tranform and pool results
logit.ci = logit.transform.pool(results, n.imp)
regular.ci = regular.pool(results, n.imp)
true.ci = rep(true.c, 3)

# Combine
study.result = rbind(logit.ci, regular.ci, true.ci)

# Forest plot
metafor::forest(slab = c('Logit', 'Regular', 'True'),
                x = study.result[,2], 
                ci.lb = study.result[,1],
                ci.ub = study.result[,3],
                refline = study.result[,2][nrow(study.result)], 
                xlab = "C-statistic")


 # IRRELEVANT ----
perf.on.validation = sapply(model.list, function(x){
    c = sapply(testlist, function(z){
        pred = predict(x, newdata=z, type = 'response')
        auc = pROC::roc(z$V3, pred, quiet = T)$auc[[1]]
        return(auc)
    })
    
    return(c)
})




## INDICATOR MATHOD TO CHECK MAR ----
testdata$missV1 = factor(ifelse(is.na(missingdata$V1), 1, 0)) # missingness conditional on V2
testdata$missV2 = factor(ifelse(is.na(missingdata$V2), 1, 0)) # missingness conditional on Y

# Fit logistic regression with indicator variable as outcome

glm(factor(missV1)~V1+V2+V3, family = 'binomial', data = testdata) %>% summary()
glm(factor(missV2)~V1+V2+V3, family = 'binomial', data = testdata) %>% summary()

testdata = testdata %>% select(V3,V2,V1)
