# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")
library(lattice)
library(mice)
library(tidyverse)

## GENERATE DATA ----
simulate.multivariate <- function(n=10000, mu=c(2, 6), sd=c(0.5, 0.8), intercept = -7){
    
    # continuous variables could be e.g., blood pressure
    continuous_var_1 = rnorm(n = n, mean = mu[1], sd = sd[1])
    continuous_var_2 = rnorm(n = n, mean = mu[2], sd = sd[2])
    # Create linear combination 
    # with or without bias as intercept?
    lin_combination = intercept + continuous_var_1 + continuous_var_2
    
    # Probability for response variable to be 1
    # Note: Due to application for logistic regression, use inverse logit function
    prob_observed_outcome = 1/(1+exp(-lin_combination))
    
    # Check that values are not approaching either 0 or 1 to avoid too deterministic approach
    summary(prob_observed_outcome)
    
    # Binary outcome variable as Bernoulli response variable
    # e.g., diabetes positive (1) or negative (0)
    # Desired probability for outcome_var = 1 between 20% and 30% to consider imbalance
    outcome_var = rbinom(n = n, size = 1, prob = prob_observed_outcome)
    
    # Combine it to a data frame
    df_complete = data.frame( 
        continuous_var_1,
        continuous_var_2,
        outcome_var
    )
    colnames(df_complete) <- c('V1', 'V2', 'V3')
    return(df_complete)

}
## AMPUTE WITH MICE ----

amputation = function(data, method='MAR', on ='covariates', prop = 0.9){
    if (method == 'MAR' && on == 'covariates'){
        ampdata = mice::ampute(data, prop = prop, patterns =  rbind(c(0,1,1), c(1,0,1)), 
                               mech = 'MAR', weights = rbind(c(0,5,0),c(0,0,5)), freq = c(.5,.5), bycases = T)$amp}
    else if (method == 'MNAR' && on == 'outcome'){
        ampdata = mice::ampute(data, prop = prop, patterns =  c(1,1,0), 
                               mech = 'MAR', bycases = T)$amp}
    else if (method == 'MAR' && on == 'outcome'){
        ampdata = mice::ampute(data, prop = prop, patterns =  c(1,1,0), 
                               mech = 'MAR', weights = rbind(c(.5,.5, 0)), bycases = T)$amp}
    else if (method == 'MNAR' && on == 'covariates'){
        ampdata = mice::ampute(data, prop = prop, patterns = rbind(c(0,1,1), c(1,0,1)), 
                               mech = 'MAR', freq=c(.5,.5), bycases = T)$amp}
    return(ampdata)
}

imputation = function(missingdata, n.imp){
    invisible(capture.output(imp <- mice(missingdata, n.imp, method = c('pmm','pmm','logreg'))))
    pred <- make.predictorMatrix(imp$data)
    pred[c("V1", "V2"), "V2"] <- 0
    invisible(capture.output(imp <- mice(imp$data, n.imp, pred = pred, maxit = n.imp+10,
                                         print = FALSE, seed = 11)))
    # add extra iterations just to be sure that the variability between the trace lines,
    # should equal the variability within each of the lines
    return(imp)
}


#define function for bootstrapping----

func = function(data, ind){
    fit = glm(factor(V3)~.,family = 'binomial', data=data[ind,])
    #pROC::auc(fit$y, fit$fitted.values, quiet = T)
    DescTools::Cstat(fit)
}

true.c.stat = function(data){
    testdata = simulate.multivariate(10000)
    fit = glm(factor(V3)~.,family = binomial, data=data)
    pred = predict(fit, newdata=testdata, type= 'response')
    auc = pROC::roc(testdata$V3, pred, quiet = T)$auc[[1]]
    return(auc)
}

# NON-PARAMETRIC BOOTSTRAP for determining within-variance of C-statistic ----
c.stat.bootstrap = function(imp, nr_boots){
    # empty matrix to store results
    # create stack of m imputed data sets
    imp_tot <- complete(imp, "long")
    # BOOTSTRAP
    # create list of imputed datasets
    imp.list <- lapply(1:imp$m, function(x){subset(imp_tot, .imp == x)})
    # apply boostrap function on every imputed data set
    bootstraps = lapply(imp.list, function(x){
        res = boot::boot(x,func, R=nr_boots)
        res$t})
    # store mean and standard deviation
    results =  cbind(sapply(bootstraps, mean), sapply(bootstraps, sd)) %>% data.frame()
    colnames(results) = c('C.statistic', 'SE')
    return(results)
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

logit.transform.pool = function(results, n.imp){
    #  Wald methode te gebruiken VAR(logitauc) =  [(logit (cub) − logit (clb)) /(2 × 1.96)]^ 2
    logit = metamisc::ccalc(cstat = results$C, cstat.cilb = results$`95%lb`,
                           cstat.ciub = results$`95%ub`, g = "log(cstat/(1-cstat))")
    # calculate logit pooled estimate
    est_c_se_log = results$SE / (results$C * (1-results$C))
    
    # wald

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

regular.pool = function(results, n.imp){
    results = metamisc::ccalc(cstat = results$C, cstat.cilb = results$`95%lb`,
                           cstat.ciub = results$`95%ub`)
    # calculate regular pooled estimate
    #VAR = ((results$ciub-results$cilb) /(2*1.96))**2
    pooled.se = pooled_se(results$theta, results$theta.se, n.imp=n.imp)
    reg.C = c('95%.Low' = mean(results$theta) - (pooled.se[2]*pooled.se[1]/sqrt(length(results))),
              Pooled.C.statistic = mean(results$theta),
              '95%.Up' = mean(results$theta) + (pooled.se[2]*pooled.se[1]/sqrt(length(results))))
    return(reg.C)
}

arcsine.pool = function(results, n.imp){
    # Pooled SE
    arcsin =  asin(sqrt(results$C.statistic))
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

# transform and pool
transform_and_pool = function(results, true.c, n.imp){
    logit = logit.transform.pool(results, n.imp)
    arcsine = arcsine.pool(results, n.imp)
    regular = regular.pool(results, n.imp)
    true = rep(true.c, 3)
    combined = rbind(regular, logit, arcsine, true) %>% data.frame()
    colnames(combined) = c('95%lb', 'C', '95%ub')
    return(combined)
}

# forestplot function
forestplot = function(study.results){
    plots = list()
    ## forest plot
    for (i in 1:length(study.results)){
        
        plots[[i]] = metafor::forest(slab = c('Logit', 'Regular', 'True'),
                                     x = study.results[[i]][,2], 
                                     ci.lb = study.results[[i]][,1],
                                     ci.ub = study.results[[i]][,3],
                                     refline = study.results[[i]][3,2], xlab = "C-statistic")
        title(paste("Simulation", i), line = -3)
    }
    return(plots)
}

# convergence rate function
convergence.rate = function(study.results){
    cr.logit = 0
    cr.regular = 0
    cr.arcsine = 0
    # does the true c lie within the intervals?
    true.C = study.results[[1]][3,2]
    for (i in 1:length(study.results)){
        
        if (true.C >= study.results[[i]][1,1] & true.C <= study.results[[i]][1,3]){
            cr.logit = cr.logit+ 1
        }
        if (true.C >= study.results[[i]][2,1] & true.C <= study.results[[i]][2,3]){
            cr.regular = cr.regular + 1
        }
    }
    # calculate proportion of times the c-statistic lies within the boostrap intervals
    cr.logit = paste0(cr.logit/length(study.results)*100, '%')
    cr.regular = paste0(cr.regular/length(study.results)*100, '%')
    return(list('regular.pooled'= cr.regular, 'logit.transformed' = cr.logit))
}

# Parallel computation simulation study
simulation.study = function(repetitions, mech, var, n.imp, nr_boots){
    #install.packages("parallel")
    library(parallel)
    # create complete data sets
    sims = replicate(repetitions, simulate.multivariate(), simplify = FALSE)
    # create clusters for parallel computing
    cl <- makeCluster(getOption("cl.cores", detectCores()))
    # load libraries and functions
    clusterEvalQ(cl, {library(mice); library(tidyverse); library(boot)})
    clusterExport(cl, c("amputation", "c.stat.bootstrap", "logit.transform.pool", 
                        "regular.pool", "pooled_se",'true.c.stat',"simulate.multivariate", "func", "arcsine.pool",
                        'n.imp','nr_boots', 'imputation'), 
                  envir=environment())
    # apply parallel computing
    study.results = parLapply(cl, sims, function(x, method=mech, on=var){
        # determine "true" c-statistic
        set.seed(11)
        true = true.c.stat(x)
        true.C = rep(true, 3)
        set.seed(11)
        # create missing data according to user-specified mechanism
        missingdata = amputation(x, method, on)
        missingdata$V3 = as.factor(missingdata$V3) # change outcome to factor
        # Impute
        set.seed(11)
        imp = imputation(missingdata, n.imp)
        # non-parametric bootstrapping
        set.seed(11)
        result = c.stat.bootstrap(imp, nr_boots = nr_boots)
        # (transform and) pool results
        logit.C = logit.transform.pool(results = result, n.imp)
        reg.C = regular.pool(results = result, n.imp)
        arcsine.C = arcsine.pool(results = result, n.imp)
        # combine pooled estimates
        combined = rbind(reg.C, logit.C, arcsine.C, true.C) %>% data.frame()
        colnames(combined) = c('95%lb', 'C', '95%ub')
        return (combined)
    })
    stopCluster(cl)
    return(study.results)
}



# Parallel computation simulation study
simulation.study.Altman = function(repetitions, mech, var, n.imp){
    #install.packages("parallel")
    library(parallel)
    # create complete data sets
    sims = replicate(repetitions, simulate.multivariate(), simplify = FALSE)
    # create clusters for parallel computing
    cl <- makeCluster(getOption("cl.cores", detectCores()))
    # load libraries and functions
    clusterEvalQ(cl, {library(mice); library(tidyverse); library(boot); library(auctestr)})
    clusterExport(cl, c("amputation", "logit.transform.pool", 
                        "regular.pool", "pooled_se",'true.c.stat',"simulate.multivariate", 
                        "arcsine.pool", 'n.imp',"transform_and_pool", 'imputation'), 
                  envir=environment())
    # apply parallel computing
    study.results = parLapply(cl, sims, function(x, method=mech, on=var){
        # determine "true" c-statistic
        # determine 'true' C-statistic first
        set.seed(11)
        true.c = true.c.stat(x)
        missingdata = amputation(x, method = 'MAR', on = 'covariates')
        missingdata$V3 = as.factor(missingdata$V3)
        # split missingdata into development and validation set
        smp = sample(seq_len(floor(nrow(missingdata)*0.5)), replace = F)
        trainset = missingdata[smp,]
        testset =  missingdata[-smp,]
        # Impute
        train.imp = imputation(trainset, n.imp)
        test.imp = imputation(testset, n.imp)
        
        # fit a model on each imputed data set (SHOULD I USE A SHRINKAGE FACTOR?)
        trainset = complete(train.imp, 'long')
        testset = complete(test.imp, 'long')
        
        # develop model by combining coefficients with Rubin's rules
        model.list = lapply(1:train.imp$m, function(x){
            d = subset(trainset, .imp == x)
            d = d %>% select(V3,V2,V1)
            fit = glm(factor(V3)~.,family = 'binomial', data=d)})
        
        # combine into a single model
        final.model = summary(pool(model.list,  rule = "rubin1987"))
        
        # Calculate predicted risk per model on complete data set
        #perf.on.complete = sapply(model.list, function(x){
        #pred = predict(x, newdata=testdata, type = 'response')
        #auc = pROC::roc(testdata$V3, pred, quiet = T)$auc[[1]]
        #return(auc)
        #})
        
        # Calculate predicted risk of final model on validation data set
        testlist <- lapply(1:test.imp$m, function(x){subset(testset, .imp == x)})
        
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
        # (Transform and) pool
        logit.ci = logit.transform.pool(results, n.imp)
        regular.ci = regular.pool(results, n.imp)
        true.ci = rep(true.c, 3)
        
        study.result = rbind(logit.ci, regular.ci, true.ci)
    })
    stopCluster(cl)
    return(study.results)
}

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
