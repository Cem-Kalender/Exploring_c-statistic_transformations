# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")
library(lattice)
library(mice)
library(tidyverse)

## GENERATE DATA ----
simulate.multivariate <- function(n=5000, mu=c(10, 6), sd=c(12, 5), intercept = 5, noise = 0){
    # Simulate
    vars = lapply(1:length(mu), function(x){rnorm(n=n,mean=mu[x],sd=sd[x])})
    vars = sapply(1:length(vars), function(x){return(vars[[x]])}) %>% as.data.frame()
    # Create linear combination
    lin_combination = sapply(1:n, function(x){rowSums(vars[x,])})
    # Probability for response variable to be 1
    # Note: Due to application for logistic regression, use inverse logit function
    prob_observed_outcome = exp(-lin_combination+intercept)/(1+exp(-lin_combination+intercept))
    # Binary outcome variable as Bernoulli response variable
    outcome_var = rbinom(n = n, size = 1, prob = prob_observed_outcome)
    # Combine it to a data frame
    vars$Y = outcome_var
    # Introduce additional randomness with noise parameter
    noise_param.1 = sample(seq_len(floor(length(outcome_var[outcome_var==1])*noise)))
    noise_param.0 = sample(seq_len(floor(length(outcome_var[outcome_var==0])*noise)))
    
    vars$Y[noise_param.0] = 1
    vars$Y[noise_param.1] = 0
    return(vars)
}
## AMPUTE WITH MICE ----
amputation = function(data){
    ampdata = mice::ampute(data, 
                           prop = 0.9, 
                           patterns =  rbind(c(0,1,1), c(1,1,0)),
                           mech = 'MAR', 
                           weights = rbind(c(0,5,0),c(5,0,0)), 
                           freq = c(.5,.5), 
                           bycases = T)$amp
    return(ampdata)
}
## Function for SE pooling ----
pooled_se <- function(est, se, n.imp){
    Qbar <- mean(est)
    U <- sum(se**2)/n.imp # within-variance
    B <- sum((est - Qbar)**2)/(n.imp-1) # between variance
    se_total <- sqrt(U + (1+1/n.imp)*B)
    r <- (1 + 1 / n.imp) * (B / U)
    v <- (n.imp - 1) * (1 + (1/r))^2
    t <- qt(0.975, v)
    res <- c(se_total, t)
    return(res)
}
# Logit transformation and pooling ----
logit.pool = function(results, n.imp){
    #  Wald methode te gebruiken VAR(logitauc) =  [(logit (cub) ??? logit (clb)) /(2 × 1.96)]^ 2
    logit = metamisc::ccalc(cstat = results$C, results$SE, g = "log(cstat/(1-cstat))")
    # calculate logit pooled estimate
    logit.pooled.se = pooled_se(logit$theta, logit$theta.se , n.imp=n.imp)
    # pooled c-statistic estimate
    c = exp(mean(logit$theta))/(1 + exp(mean(logit$theta)))
    # pooled SE
    se = exp(logit.pooled.se[1])/(1+exp(logit.pooled.se[1]))
    return(c(c, se, logit.pooled.se[2]))
}
# Regular pooling ----
regular.pool = function(results, n.imp){
    pooled.se = pooled_se(results$C, results$SE, n.imp=n.imp)
    reg.C = c('95%.Low' = mean(results$C) - (pooled.se[2]*(pooled.se[1]/sqrt(length(results)))),
              Pooled.C.statistic = mean(results$C),
              '95%.Up' = mean(results$C) + (pooled.se[2]*(pooled.se[1]/sqrt(length(results)))))
    return(c(mean(results$C), pooled.se[1], pooled.se[2]))
}
# Arcsine transformation and pooling ----
arcsine.pool = function(results, n.imp){
    # Pooled SE
    arcsin =  asin(sqrt(results$C))
    est_c_se_arc =asin(sqrt(results$SE))
    arcsin.pooled.se = pooled_se(arcsin, est_c_se_arc, n.imp=n.imp)
    # return pooled c-statistic estimate and pooled se
    arc = c(sin(mean(arcsin))**2, sin(arcsin.pooled.se[1])**2)
    return(c(arc, arcsin.pooled.se[2]))
}
# SE AUC (wilcoxon)
se_auc <- function(auc, n_p, n_n) {
    D_p = (n_p - 1) * ((auc/(2 - auc)) - auc^2)
    D_n = (n_n - 1) * ((2 * auc^2)/(1 + auc) - auc^2)
    SE_auc = sqrt((auc * (1 - auc) + D_p + D_n)/(n_p * n_n))
    return(SE_auc)
}
# Coverage process
coverage.process = function(study.results){
    lapply(1:3, function(x){
        lapply(1:length(study.results[[1]]), function(j){
            # logit
            z = study.results[[x]][[j]][1,3]
            logit = log(study.results[[x]][[j]][1,c(1,2)]/(1-study.results[[x]][[j]][1,c(1,2)]))
            # lower bound
            logit.ub = exp(logit[1] + z*logit[2]) /
                (1 + exp(logit[1] + z*logit[2]))
            # upper bound
            logit.lb =  exp(logit[1] - z*logit[2]) /
                (1 + exp(logit[1] - z*logit[2]))
            # combine
            logit.C = c('95%.Low' = logit.lb,
                        Pooled.C.statistic = exp(logit[1])/(1 + exp(logit[1])),
                        '95%.Up' =logit.ub)
            # regular
            z = study.results[[x]][[j]][2,3]
            reg.C = c('95%.Low' = study.results[[x]][[j]][2,1] - (z*study.results[[x]][[j]][2,2]),
                      Pooled.C.statistic = study.results[[x]][[j]][2,1],
                      '95%.Up' = study.results[[x]][[j]][2,1] + (z*study.results[[x]][[j]][2,2]))
            
            # arcsine
            z = study.results[[x]][[j]][3,3]
            arcsin = asin(sqrt(study.results[[x]][[j]][3,c(1,2)]))
            # lower bound
            arcsin.lb = sin(arcsin[1] - (z*arcsin[2]))**2
            # Upper bound
            arcsin.ub = sin(arcsin[1] + (z*arcsin[2]))**2
            
            # combine
            arcsin.C = c('95%.Low' = arcsin.lb,
                         Pooled.C.statistic = sin(arcsin[1])**2,
                         '95%.Up' = arcsin.ub)
            # combine all
            CI = rbind(logit.C, reg.C, arcsin.C,'true.c'=rep(study.results[[x]][[j]][4,1], 3))
            
            return(CI)
        })
    })
}
# Convergence rate function ----
coverage.rate = function(study.results){
   
    convergence = lapply(1:3, function(i){
        cr.logit = 0
        cr.regular = 0
        cr.arcsine = 0
        
        convergence = sapply(1:length(study.results[[1]]), function(j){
            # does the true c lie within the intervals?
            true.C = study.results[[i]][[j]][4,2]
            if (true.C >= study.results[[i]][[j]][1,1] & true.C <= study.results[[i]][[j]][1,3]){
                cr.logit = cr.logit + 1
            }
            if (true.C >= study.results[[i]][[j]][2,1] & true.C <= study.results[[i]][[j]][2,3]){
                cr.regular = cr.regular + 1
            }
            if (true.C >= study.results[[i]][[j]][3,1] & true.C <= study.results[[i]][[j]][3,3]){
                cr.arcsine = cr.arcsine + 1}
            
            c('regular'=cr.regular, 'logit'=cr.logit, 'arcsine'=cr.arcsine)
        })
    })
    convergence = sapply(convergence, function(i){rowSums(i)})/length(study.results[[1]])
    
    return(convergence)
}
# Bias ----
bias = function(results){
    logit.bias = sapply(results, function(x){
        c = x[4,1]
        estimate = x[1,1]
        logit.bias = estimate - c
        return(logit.bias)
    })
    
    regular.bias = sapply(results, function(x){
        c = x[4,1]
        estimate = x[2,1]
        regular.bias = estimate - c
        return(regular.bias)})
    
    arcsine.bias = sapply(results, function(x){
        c = x[4,1]
        estimate = x[3,1]
        arcsine.bias = estimate - c
        return(arcsine.bias)})
    
    bias = cbind(regular.bias = sum(unlist(regular.bias)), 
                 logit.bias = sum(unlist(logit.bias)), 
                 arcsine.bias = sum(unlist(arcsine.bias)))/length(results)
    return(bias)
}
# Parallel computation simulation study (non) ----
simulation.study = function(repetitions, imp.method, n.imp, noise_param){
    # Install.packages("parallel")
    library(parallel)
    # Create complete data sets
    sims = replicate(repetitions, simulate.multivariate(noise=noise_param), simplify = FALSE)
    # Create clusters for parallel computing
    cl <- makeCluster(getOption("cl.cores", detectCores()))
    # Load libraries and functions
    clusterEvalQ(cl, {library(mice); library(tidyverse); library(boot)})
    clusterExport(cl, c("amputation",'logit.pool', "regular.pool",
                        "pooled_se","simulate.multivariate", "arcsine.pool", 'n.imp', 'se_auc'), 
                  envir=environment())
    # Apply parallel computing
    study.results = parLapply(cl, sims, function(x, imp_method=imp.method,noise=noise_param){
        ################### DETERMINE THETA ############
        sm = sample(seq_len(floor(nrow(x)*0.5)), replace = F)
        testdata = x[sm,] # TRUE or FALSE?
        data = x[-sm,]
        fit = glm(factor(Y)~.,family = 'binomial', data=data)
        pred = predict(fit, newdata=testdata, type= 'response')
        true.c = pROC::roc(testdata$Y, pred, quiet = T)$auc[[1]]
        
        ###################### AMPUTE ########################
        missingdata = amputation(x)
        missingdata$Y = as.factor(missingdata$Y)
        
        # Split
        trainset =  missingdata[-sm,]
        testset = missingdata[sm,]
        ####################### IMPUTE #######################
      
        train.imp = mice(trainset, method = c('norm','norm','logreg'), n.imp)
        test.imp = mice(testset, method = c('norm','norm','logreg'), n.imp)
        
        # fit a model on each imputed data set
        trainset = complete(train.imp, 'long')
        testset = complete(test.imp, 'long')
        
        # Develop model by combining coefficients with Rubin's rules
        model.list = lapply(1:train.imp$m, function(x){
            d = subset(trainset, .imp == x)
            d = d %>% select(Y,V2,V1)
            fit = glm(factor(Y)~.,family = 'binomial', data=d) # Shrinkage?
        })
        ########################### POOL MODELS #######################
        final.model = summary(mice::pool(model.list,  rule = "rubin1987"))
        
        # Make subsets of validation data set
        testlist <- lapply(1:test.imp$m, function(x){subset(testset, .imp == x)})
        
        # Prepare data
        validation.data <- lapply(testlist, function(x){
            model.matrix(formula(paste0("~", c('V2+V1'))), data=x)})
        # Predict
        predictions = lapply(validation.data, function(x){
            pred = sapply(1:nrow(x), function(i){
                1/(1+exp(-x[i,] %*%  final.model[c(1,2,3),2]))
            })
        })
        ################## COMPUTE THETA HAT #################
        perf.on.validation = lapply(1:n.imp, function(x){
            pred = predictions[[x]]
            test = testlist[[x]]
            auc = pROC::roc(test$Y, pred, quiet = T)$auc[[1]]
            se.auc = se_auc(auc, n_p = sum(test$Y==1), n_n = sum(test$Y==0))
            auc_and_se = cbind('C'=auc, 'SE'=se.auc)
            return(auc_and_se)
        })
        
        # Results
        results = data.frame(matrix(unlist(perf.on.validation), nrow=length(perf.on.validation), byrow=TRUE))
        colnames(results) <- c('C', 'SE')
        
        ################## TRANSFORM AND POOL ###################
        logit.ci = logit.pool(results, n.imp)
        regular.ci = regular.pool(results, n.imp)
        arcsine.ci = arcsine.pool(results, n.imp)
        true.ci = rep(true.c, 2)
        
        # Combine
        study.result = rbind('logit'=logit.ci, 'reg'=regular.ci, 'arc'=arcsine.ci, true.ci)
        
    })
    stopCluster(cl)
    return(study.results)
}
# convenient function to turn lists into data frames ----
df.maker = function(l){
    return(data.frame(matrix(unlist(l), nrow=length(l), byrow=TRUE)))
}