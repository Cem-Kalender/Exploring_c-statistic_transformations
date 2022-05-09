# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")
library(lattice)
library(mice)
library(tidyverse)

## GENERATE DATA ----
simulate.multivariate <- function(n=10000, mu=c(22, 12, 0), threshold = NULL, sigma = c(1.0, 0.3, 0.5, 
                                                                                        0.3, 1.0, 0.5, 
                                                                                        0.5, 0.5, 1.0), prob=0.5, nr_conditions = 2){
    sigma = matrix(sigma,3,3)
    data = MASS::mvrnorm(n = n,mu = mu,Sigma = sigma) %>% as.data.frame()
    
    # Define conditions to create a binary variable
    if (nr_conditions==1){
        
        condition = data$V1>median(data$V1) & data$V2>median(data$V2)
        # generate binomial variable for V3
        data$V3 = ifelse(condition,rbinom(nrow(data[condition,]), size = 1, prob = prob), 0)}
    
    else if (nr_conditions==2){
        
        condition = data$V1>median(data$V1) & data$V2>median(data$V2)
        condition2 = data$V1<median(data$V1) & data$V2>median(data$V2)
        # generate binomial variable for V3
        data$V3 = ifelse(condition,rbinom(nrow(data[condition,]), size = 1, prob = prob),
                         ifelse(condition2,rbinom(nrow(data[condition2,]), size = 1, prob = 0.1), 0))}
    
    else if (nr_conditions == 3){
        
        condition = data$V1>median(data$V1) & data$V2>median(data$V2)
        condition2 = data$V1<median(data$V1) & data$V2>median(data$V2)
        condition3 = data$V1>median(data$V1) & data$V2<median(data$V2)
        # generate binomial variable for V3
        data$V3 = ifelse(condition,rbinom(nrow(data[condition,]), size = 1, prob = prob),
                         ifelse(condition2,rbinom(nrow(data[condition2,]), size = 1, prob = 0.1),
                                ifelse(condition3,rbinom(nrow(data[condition3,]), size = 1, prob = 0.1),0)))
    }
    
    else if (nr_conditions == 'linear'){
        # linear model with sigmoid transformation of linear predictors
        mod = lm(V3~V2+V1, data = data)
        y = 1/(1+exp(mod$fitted.values))
        # introduce randomness/noise by means of bernoulli trial
        data$V3 = ifelse(y>threshold,rbinom(nrow(data[y>threshold,]), size = 1, prob), 0)
        data$V3[data$V3==0] = rbinom(length(data$V3[data$V3==0]), size = 1, prob = 1-prob)}
    
    else {
        stop('Conditions not correctly specified')
    }
    return(data)
}

## AMPUTE WITH MICE ----

amputation = function(data, pattern_list, prop_list, mech_list, probabilities){
    # if the length of the inputs are unequal, return NA
    len_inputs = sapply(list(pattern_list, prop_list, mech_list, probabilities), length) %>% unique()
    if (length(len_inputs) == 1){
        # empty list to store amputations
        amputations = list()
        for (i in 1:length(pattern_list)){
            amputations[[i]]=mice::ampute(data, patterns = pattern_list[[i]], prop = prop_list[i], mech = mech_list[i])$amp
        }
        # sample from amputations with probability
        indices <- sample(x = seq_len(length(pattern_list)), size = nrow(data), 
                          replace = TRUE, prob = probabilities)
        # empty matrix to store
        ampdata <- matrix(NA, nrow = nrow(data), ncol = ncol(data)) %>% as.data.frame()
        # insert amputations based on index
        for (i in 1:length(pattern_list)){
            ampdata[indices == i, ] <- as.matrix(amputations[[i]][indices == i, ])}
    }
    else {
        stop('The lengths of the inputs are not equal')
    }
    return(ampdata)
}

#define function for bootstrapping
func = function(data,ind){
    fit = glm(V3~V1+V2,family = 'binomial', data=data[ind,])
    pROC::auc(pROC::roc(fit$y, fit$fitted.values, quiet = T))
    #DescTools::Cstat(fit)
}

# NON-PARAMETRIC BOOTSTRAP for determining within-variance of C-statistic ----
C_stat.bootstrap = function(imp, Y, nr_boots){
    # empty matrix to store results
    results = matrix(nrow = imp$m, ncol = 2) %>% data.frame()
    colnames(results) = c('C.statistic', 'SE')
    # create stack of m imputed data sets
    imp_tot <- complete(imp, "stacked")
    start = seq(1, n.imp*nrow(missingdata), by=nrow(missingdata))
    end = seq(nrow(missingdata), n.imp*nrow(missingdata), by=nrow(missingdata))
    # BOOTSTRAP
    for (i in 1:imp$m){
        # data set m
        data.copy = imp_tot[start[i]:end[i],]
        ## create bootstrap
        res = boot::boot(data.copy, func, R=nr_boots)
        print(paste('Iteration', i, 'of', n.imp, 'completed'))
        # store sample statistic and corresponding error (SD=SE)
        results[i, 1] = res$t0
        results[i, 2] = sd(res$t)
    }
    return(results)
}
## Function for SE pooling ----
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