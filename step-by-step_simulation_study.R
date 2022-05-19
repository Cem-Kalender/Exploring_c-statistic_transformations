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


## INDICATOR MATHOD TO CHECK MAR ----
testdata$missV1 = factor(ifelse(is.na(missingdata$V1), 1, 0)) # missingness conditional on V2
testdata$missV2 = factor(ifelse(is.na(missingdata$V2), 1, 0)) # missingness conditional on Y

# Fit logistic regression with indicator variable as outcome

glm(factor(missV1)~V1+V2+V3, family = 'binomial', data = testdata) %>% summary()
glm(factor(missV2)~V1+V2+V3, family = 'binomial', data = testdata) %>% summary()

testdata = testdata %>% select(V3,V2,V1)
