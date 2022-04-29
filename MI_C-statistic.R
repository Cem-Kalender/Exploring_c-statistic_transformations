# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")

# import packages
library(lattice)
library(mice)
library(tidyverse)

## simulate dataset----

set.seed(2015)
#options(scipen=999)

# randomly generate a dataset to be our complete dataset


testdata <- as.data.frame(MASS::mvrnorm(n = 10000, 
                                        mu = c(22, 12, 0),      # V1   V2   V3
                                        Sigma = matrix(data = c(1.0, 0.6, 0.4, 
                                                                0.6, 1.0, 0.5, 
                                                                0.4, 0.5, 1.0), 
                                                       nrow = 3, 
                                                       byrow = T)))
summary(testdata)

# People falling in a certain region have a a higher probabilityof Y = 1
condition = testdata$V1>median(testdata$V1) & testdata$V2>median(testdata$V2)
condition2 = testdata$V1<median(testdata$V1) & testdata$V2>median(testdata$V2)
# plot
plot(testdata$V1, testdata$V2, col = ifelse(condition, 'red', 'green'))

# generate binomial variable for V3
testdata$V3 = ifelse(condition,rbinom(nrow(testdata[condition,]), size = 1, prob = 0.5),
                     ifelse(condition2,rbinom(nrow(testdata[condition2,]), size = 1, prob = 0.1), 0))
plot(testdata$V1, testdata$V2, col = ifelse(testdata$V3==1,'red', 'green'))


# visualize missing mechanism plotting data points----
par(mfrow=c(1,4))
# no missingness
plot(testdata$V1, testdata$V2, main = 'Normal')
abline(lm(testdata$V2~testdata$V1), col = 'blue')
# MAR
plot(testdata$V1, testdata$V2, col = ifelse(testdata$V1>median(testdata$V1),'red','black'), main = 'Missingness: MAR')
abline(lm(testdata$V2[testdata$V1<median(testdata$V1)]~testdata$V1[testdata$V1<median(testdata$V1)], data = testdata), col = "blue")
# MNAR
plot(testdata$V1, testdata$V2, col = ifelse(testdata$V2>median(testdata$V2),'red','black'), main = 'Missingness: MNAR')
abline(lm(testdata$V2[testdata$V2<median(testdata$V2)]~testdata$V1[testdata$V2<median(testdata$V2)], data = testdata), col = "blue")

# create MAR+MNAR missingness mechanism and ampute ----

# first split the complete data, before amputing
samp = sample(seq_len(floor(nrow(testdata)*0.8)))
test.data = testdata[-samp,] # this the test data that has not been amputed

# we will now ampute
missingdata = testdata[samp,]

# defining missing mechanisms
missing.mechanism1 = missingdata$V2<median(missingdata$V2) & missingdata$V1<median(missingdata$V1)
missing.mechanism2 = missingdata$V2>median(missingdata$V2) & missingdata$V1>median(missingdata$V1)
missinglist = ifelse(missing.mechanism1, rbinom(nrow(missingdata[missing.mechanism1,]), 1, 0.4),
                     ifelse(missing.mechanism2, rbinom(nrow(missingdata[missing.mechanism1,]), 1, 0.7),0))

# visualize
plot(missingdata$V1, missingdata$V2, col = ifelse(missinglist==1, 'red', 'green'),
     main = 'Incorporated missing mechanism')
abline(lm(testdata$V2~testdata$V1), col = "blue")

# Ampute (univariate missingness)
missingdata$V1[missinglist==1] = NA
abline(lm(missingdata$V2~missingdata$V1), col = "red")

sum(is.na(missingdata$V1)) # 3100 missing values created

## MI ----

#data = result$amp

N = sum(is.na(missingdata$V1))
miss = data[is.na(missingdata$V1),]
n.imp = 100
results = matrix(nrow = n.imp, ncol = 3) %>% data.frame()
colnames(results) = c('C.statistic', '95%CI.lower', '95%CI.upper')


for (i in 1:n.imp){
    # copy missing data
    data.copy = missingdata
    
    # create imputation model
    imp.model = lm(V1~factor(V3, c(0,1))+V2, data=data.copy)
    
    #drawing coefficients from a multivariate Student-t distribution
    betas <- rmt(N, mean=imp.model$coef, S=vcov(imp.model), df=imp.model$df.residual)
    
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
    
    # Use imputed data to fit a prediction model
    pred.model = glm(factor(V3, c(0,1))~V2+V1, family = 'binomial', data.copy)
    # predict
    preds = predict(pred.model, type = 'response', newdata = test.data)
    
    # calculate C-statistic and corresponding 95%CI 
    ROC = pROC::roc(test.data$V3, preds)
    C <- pROC::auc(ROC)
    ci = pROC::ci(C)
    
    # store estimates and CI
    results[i, 1] <-  C
    results[i, c(2,3)] <- ci
    
}

# inspect results
summary(results)

par(mfrow=c(1,3))
# plot histogram
hist(results$C.statistic, freq = F, breaks = 20, main = 'Distribution of C-statistic')
# plot densityplot
dens <- density(results$C.statistic)
lines(dens, col = 'red')

#logit transformation of the c-statistics
logit.C =  log(results$C.statistic/(1 - results$C.statistic))

# plot histogram
hist(logit.C, freq = F, breaks = 20, main = 'Distribution of logit transformed C-statistic')
# plot densityplot
dens <- density(logit.C)
lines(dens, col = 'red')

#logit transformation of the c-statistics
arcsine.C =  asin(sqrt(results$C.statistic))

# plot histogram
hist(arcsine.C, freq = F, breaks = 20, main = 'Distribution of arcsine transformed C-statistic')
# plot densityplot
arcsine.dens <- density(arcsine.C)
lines(arcsine.dens, col = 'red')










## Rubin's rules----


## pooled estimates is just the average:
beta = results$Estimate %>% mean()

# The SE for the pooled regression coefficient estimate is retrieved using Rubin's rules
Qbar <- results$Estimate %>% mean()
U <- sum(results$SE**2)/n.imp
B <- sum((results$Estimate - Qbar)**2)/(n.imp-1)
se.beta <- sqrt(U + (1+1/n.imp)*B)
# Inspect
c(beta=beta,se.beta=se.beta)


# univariate amputation with 10 percent missing----
result <- ampute(testdata, prop = 0.1, patterns=c(0,1,1), mech = 'MAR')
result$amp$V1
# inspect missingness patterns
mypatterns <- result$patterns
mypatterns


# multivariate amputation----
# add missingesness patterns
mypatterns[2, 1] <- 0
mypatterns <- rbind(mypatterns, c(0, 1, 0))

# define frequency of missingness patterns
myfreq <- c(0.7, 0.1, 0.1, 0.1)

# define weights
myweights <- result$weights
myweights[1, ] <- c(0, 0.8, 0.4)
myweights[3, ] <- c(3.0, 1.0, 0)

#the type of logistic probability distribution that is applied to the weighted sum scores

# Ampute again

result <- ampute(testdata, freq = myfreq, 
                 patterns = mypatterns, mech = "MAR")
result

