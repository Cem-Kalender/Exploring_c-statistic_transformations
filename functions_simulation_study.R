# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")
library(lattice)
library(mice)
library(tidyverse)

## GENERATE DATA ----
simulate.multivariate <- function(n=5000, mu=c(10, 6), sd=c(12, 5), intercept = -5, noise = 0){
    # Simulate
    vars = lapply(1:length(mu), function(x){rnorm(n=n,mean=mu[x],sd=sd[x])})
    vars = sapply(1:length(vars), function(x){return(vars[[x]])}) %>% as.data.frame()
    # Create linear combination
    lin_combination = sapply(1:n, function(x){intercept + rowSums(vars[x,])})
    # Probability for response variable to be 1
    # Note: Due to application for logistic regression, use inverse logit function
    prob_observed_outcome = 1/(1+exp(lin_combination))
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
amputation = function(data, method='MAR', on ='X&Y', prop = 0.9){
    if (method == 'MAR' && on == 'X&Y'){
        ampdata = mice::ampute(data, prop = prop, patterns =  rbind(c(0,1,1), c(1,1,0)), 
                               mech = 'MAR', weights = rbind(c(0,5,0),c(5,0,0)), freq = c(.5,.5), bycases = T)$amp}
    else if (method == 'MNAR' && on == 'outcome'){
        ampdata = mice::ampute(data, prop = prop, patterns =  c(1,1,0), 
                               mech = 'MAR', bycases = T)$amp}
    else if (method == 'MAR' && on == 'outcome'){
        ampdata = mice::ampute(data, prop = prop, patterns =  c(1,1,0), 
                               mech = 'MAR', weights = rbind(c(5,5,0)), bycases = T)$amp}
    else if (method == 'MNAR' && on == 'covariates'){
        ampdata = mice::ampute(data, prop = prop, patterns = rbind(c(0,1,1), c(1,0,1)), 
                               mech = 'MAR', freq=c(.5,.5), bycases = T)$amp}
    return(ampdata)
}
## IMPUTE WITH MICE ----
imputation = function(missingdata, imp_method, n.imp){
    invisible(capture.output(imp <- mice(missingdata, maxit = 0)))
    pred <- imp$pred
    pred[c('V1', 'V2'),'V2']=0
    pred[c('V1', 'V2'),'V1']=0
    
    invisible(capture.output(imp <- mice(imp$data, n.imp, pred = pred, method=imp_method,  maxit = n.imp,
                                         print = FALSE)))
    # add extra iterations just to be sure that the variability between the trace lines,
    # should equal the variability within each of the lines
    return(imp)
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
logit.transform.pool = function(results, n.imp){
    #  Wald methode te gebruiken VAR(logitauc) =  [(logit (cub) ??? logit (clb)) /(2 × 1.96)]^ 2
    logit = metamisc::ccalc(cstat = results$C, results$SE, g = "log(cstat/(1-cstat))")
    # calculate logit pooled estimate
    est_c_se_log = results$SE / (results$C * (1-results$C))
    logit.pooled.se = pooled_se(logit$theta, est_c_se_log , n.imp=n.imp)
    
    logit.ub = exp(mean(logit$theta) + (logit.pooled.se[2]*logit.pooled.se[1])) /
        (1 + exp(mean(logit$theta) + (logit.pooled.se[2]*logit.pooled.se[1])))
    
    logit.lb = exp(mean(logit$theta) - (logit.pooled.se[2]*logit.pooled.se[1])) /
        (1 + exp(mean(logit$theta) - (logit.pooled.se[2]*logit.pooled.se[1])))
    
    logit.C = c('95%.Low' = logit.lb,
                Pooled.C.statistic = exp(mean(logit$theta))/(1 + exp(mean(logit$theta))),
                '95%.Up' =logit.ub)
    return(logit.C)
}
# Regular pooling ----
regular.pool = function(results, n.imp){
    results = metamisc::ccalc(cstat = results$C, results$SE)
    # calculate regular pooled estimate
    pooled.se = pooled_se(results$theta, results$theta.se, n.imp=n.imp)
    reg.C = c('95%.Low' = mean(results$theta) - (pooled.se[2]*pooled.se[1]/sqrt(length(results))),
              Pooled.C.statistic = mean(results$theta),
              '95%.Up' = mean(results$theta) + (pooled.se[2]*pooled.se[1]/sqrt(length(results))))
    return(reg.C)
}
# Arcsine transformation and pooling ----
arcsine.pool = function(results, n.imp){
    # Pooled SE
    arcsin =  asin(sqrt(results$C))
    est_c_se_arc =asin(sqrt(results$SE))
    
    arcsin.pooled.se = pooled_se(arcsin, est_c_se_arc, n.imp=n.imp)
    # Upper bound
    arcsin.ub = sin(mean(arcsin))**2 + sin((arcsin.pooled.se[2]*arcsin.pooled.se[1]/sqrt(length(results))))**2
    arcsin.ub = ifelse(arcsin.ub>1, 1, arcsin.ub)
    # lower bound
    arcsin.lb = sin(mean(arcsin))**2 - sin((arcsin.pooled.se[2]*arcsin.pooled.se[1]/sqrt(length(results))))**2
    arcsin.lb = ifelse(arcsin.lb<0, 0, arcsin.lb)
    # Inspect
    arcsin.C = c('95%.Low' = arcsin.lb,
                 Pooled.C.statistic = sin(mean(arcsin))**2,
                 '95%.Up' = arcsin.ub)
    return(arcsin.C)
}
# Forestplot function ----
forestplot = function(study.results){
    plots = list()
    ## forest plot
    for (i in 1:length(study.results)){
        
        plots[[i]] = metafor::forest(slab = c('Logit', 'Regular', 'Arcsine', 'True'),
                                     x = study.results[[i]][,2], 
                                     ci.lb = study.results[[i]][,1],
                                     ci.ub = study.results[[i]][,3],
                                     refline = study.results[[i]][4,2], xlab = "C-statistic")
        title(paste("Simulation", i), line = -3)
    }
    return(plots)
}
# Convergence rate function ----
coverage.rate = function(study.results){
    cr.logit = 0
    cr.regular = 0
    cr.arcsine = 0
    convergence = sapply(1:length(study.results), function(i){
        
        # does the true c lie within the intervals?
        true.C = study.results[[i]][4,2]
        if (true.C >= study.results[[i]][1,1] & true.C <= study.results[[i]][1,3]){
            cr.logit = cr.logit + 1
        }
        if (true.C >= study.results[[i]][2,1] & true.C <= study.results[[i]][2,3]){
            cr.regular = cr.regular + 1
        }
        if (true.C >= study.results[[i]][3,1] & true.C <= study.results[[i]][3,3]){
            cr.arcsine = cr.arcsine + 1}
        
        return(c('regular'=cr.regular, 'logit'=cr.logit, 'arcsine'=cr.arcsine))
    })
    convergence = data.frame(t(convergence))
    convergence = colSums(convergence)/length(study.results)
    
    return(convergence)
}
# Bias ----
bias = function(results){
    logit.bias = sapply(results, function(x){
        c = x[4,2]
        estimate =x[1,2]
        logit.bias = estimate - c
        return(logit.bias)
    })
    
    regular.bias = sapply(results, function(x){
        c = x[4,2]
        estimate = x[2,2]
        regular.bias = estimate - c
        return(regular.bias)})
    
    arcsine.bias = sapply(results, function(x){
        c = x[4,2]
        estimate = x[3,2]
        arcsine.bias = estimate - c
        return(arcsine.bias)})
    bias = cbind(regular.bias = regular.bias, 
                 logit.bias = logit.bias, 
                 arcsine.bias = arcsine.bias)/length(results)
    return(bias)
}
# Parallel computation simulation study (non) ----
simulation.study = function(repetitions, mech, var, imp.method, n.imp, noise_param){
    # Install.packages("parallel")
    library(parallel)
    # Create complete data sets
    sims = replicate(repetitions, simulate.multivariate(noise=noise_param), simplify = FALSE)
    # Create clusters for parallel computing
    cl <- makeCluster(getOption("cl.cores", detectCores()))
    # Load libraries and functions
    clusterEvalQ(cl, {library(mice); library(tidyverse); library(boot); library(auctestr)})
    clusterExport(cl, c("amputation", "logit.transform.pool",'arcsine.pool', "regular.pool",
                        "pooled_se","simulate.multivariate", "test.performance",
                        "arcsine.pool", 'n.imp', 'imputation'), 
                  envir=environment())
    # Apply parallel computing
    study.results = parLapply(cl, sims, function(x, method=mech, on=var, imp_method=imp.method,noise=noise_param){
        
        # determine 'true' C-statistic first
        sm = sample(seq_len(floor(nrow(x)*0.5)), replace = F)
        testdata = x[sm,] # TRUE or FALSE?
        data = x[-sm,]
        fit = glm(factor(Y)~.,family = 'binomial', data=data)
        pred = predict(fit, newdata=testdata, type= 'response')
        true.c = pROC::roc(testdata$Y, pred, quiet = T)$auc[[1]]
        
        ###################### AMPUTE ########################
        trainset = amputation(data, method = 'MAR', on = 'X&Y')
        trainset$Y = as.factor(trainset$Y)
        testset =  amputation(testdata, method = 'MAR', on = 'X&Y')
        testset$Y = factor(testset$Y)
        # Split
        
        ####################### IMPUTE #######################
        train.imp = imputation(trainset, imp_method = c('norm','norm','logreg'), n.imp)
        test.imp = imputation(testset, imp_method = c('norm','norm','logreg'), n.imp)
        
        # fit a model on each imputed data set
        trainset = complete(train.imp, 'long')
        testset = complete(test.imp, 'long')
        
        # Develop model by combining coefficients with Rubin's rules
        model.list = lapply(1:train.imp$m, function(x){
            d = subset(trainset, .imp == x)
            d = d %>% select(Y,V2,V1)
            fit = glm(factor(Y)~.,family = 'binomial', data=d) # Shrinkage?
        })
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
        # Compute perfomance
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
        
        # Tranform and pool results
        logit.ci = logit.transform.pool(results, n.imp)
        regular.ci = regular.pool(results, n.imp)
        arcsine.ci = arcsine.pool(results, n.imp)
        true.ci = rep(true.c, 3)
        
        # Combine
        study.result = rbind(logit.ci, regular.ci, arcsine.ci, true.ci)
        
    })
    stopCluster(cl)
    return(study.results)
}