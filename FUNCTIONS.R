# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")
library(lattice)
library(mice)
library(tidyverse)

## GENERATE DATA ----
simulate.multivariate <- function(n=1000, mu=c(22, 12, 0), threshold = 0.6, sigma = c(1.0, 0.3, 0.5, 
                                                                                        0.3, 1.0, 0.5, 
                                                                                        0.5, 0.5, 1.0), prob=0.9){
    #install.packages("MASS")
    sigma = matrix(sigma,3,3)
    data = MASS::mvrnorm(n = n,mu = mu,Sigma = sigma) %>% as.data.frame()
    # introduce noise to handle multicollinearity
    
    # linear model with sigmoid transformation of linear predictors
    mod = lm(V3~V2+V1, data = data)
    y = 1/(1+exp(mod$fitted.values))
    # introduce randomness/noise by means of bernoulli trial
    data$V3 = ifelse(y>threshold,rbinom(nrow(data[y>threshold,]), size = 1, prob), 0)
    data$V3[data$V3==0] = rbinom(length(data$V3[data$V3==0]), size = 1, prob = 1-prob)
    # binary outcome variable as as factor
    return(data)
}
## AMPUTE WITH MICE ----

amputation = function(data, method='MAR', on ='covariates', prop = 0.9){
    if (method == 'MAR' && on == 'covariates'){
        ampdata = mice::ampute(data, prop = prop, patterns =  rbind(c(0,1,1), c(1,0,1)), 
               mech = 'MAR', weights = rbind(c(0,1,0),c(0,0,1)), freq = c(.5,.5), bycases = T)$amp}
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

#define function for bootstrapping----

func = function(data, ind){
    fit = glm(V3~V1+V2,family = 'binomial', data=data[ind,])
    #pROC::auc(fit$y, fit$fitted.values, quiet = T)
    DescTools::Cstat(fit)
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
    logit =  log(results/(1 - results))
    # calculate logit pooled estimate
    est_c_se_log = results$SE / (results$C.statistic * (1-results$C.statistic))
    logit.pooled.se = pooled_se(logit$C.statistic, est_c_se_log, n.imp=n.imp)
    logit.ub = exp(mean(logit$C.statistic) + (logit.pooled.se[2]*logit.pooled.se[1])) /
        (1 + exp(mean(logit$C.statistic) + (logit.pooled.se[2]*logit.pooled.se[1])))
    
    logit.lb = exp(mean(logit$C.statistic) - (logit.pooled.se[2]*logit.pooled.se[1])) /
        (1 + exp(mean(logit$C.statistic) - (logit.pooled.se[2]*logit.pooled.se[1])))
    
    logit.C = c('95%.Low' = logit.lb,
                Pooled.C.statistic = exp(mean(logit$C.statistic))/(1 + exp(mean(logit$C.statistic))),
                '95%.Up' =logit.ub)
    return(logit.C)
}

regular.pool = function(results, n.imp){
    # calculate regular pooled estimate
    pooled.se = pooled_se(results$C.statistic, results$SE, n.imp=n.imp)
    reg.C = c('95%.Low' = mean(results$C.statistic) - (pooled.se[2]*pooled.se[1]/sqrt(length(results))),
              Pooled.C.statistic = mean(results$C.statistic),
              '95%.Up' = mean(results$C.statistic) + (pooled.se[2]*pooled.se[1]/sqrt(length(results))))
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

# simulation study function (50% multivariate MAR/ imputation pmm)
simulation.study.loop = function(repetitions, n.imp, imp.method='pmm', nr_boots, mech = 'MAR', var = 'covariates'){
    #DEFAULT: Multivariate MAR
    
    # simulate datasets
    set.seed(11)
    sims = replicate(repetitions, simulate.multivariate(), simplify = FALSE)
    # empy list to store results
    study.results = list()
    # perform study in every simulated data set
    for (i in 1:length(sims)){
        
        cat(crayon::blue(paste('SIMULATION',i,'\n\n')))
        
        # Ampute
        cat(crayon::red('Amputing...\n'))
        set.seed(11)
        missingdata = amputation(sims[[i]], method=mech, on=var)
        missingdata$V3 = as.factor(missingdata$V3) # change outcome to factor
        
        # impute with predictive mean matching and logistic regression
        cat(crayon::red('Imputing...\n\n'))
        set.seed(11)
        invisible(capture.output(imp <- mice::mice(missingdata, n.imp, method = c(imp.method, imp.method, 'logreg'), seed = 5)))
     
        # calculate c stat complete data
        cat(crayon::white('Calculating true c-statistic...\n'))
        set.seed(11)
        true = boot::boot(sims[[i]], func, R=100)
        true.C.CI = boot::boot.ci(true, type = 'perc')$percent[c(4,5)]
        true.C = c(true.C.CI[1], true$t0, true.C.CI[2])
        
        # calculate bootstrapped c stat
        cat(crayon::white('Calculating bootstrapped c-statistic...\n\n'))
        set.seed(11)
        results = c.stat.bootstrap(imp, nr_boots = nr_boots)
    
        # regular c-stat
        logit.C = logit.transform.pool(results)
        reg.C = regular.pool(results)
        
        # combine estimates
        cat(crayon::green('Combining results...\n\n'))
        combined = rbind(logit.C, reg.C, true.C) %>% data.frame()
        colnames(combined) = c('95%lb', 'C', '95%ub')
        study.results[[i]] = combined
    }
    return (study.results)
}

# forestplot function
forestplot = function(study.results){
    plots = list()
    ## forest plot
    for (i in 1:length(study.results)){

        plots[[i]] = metafor::forest(slab = c('Regular','Logit', 'Arcsine', 'True'),
                                     x = study.results[[i]]$C, 
                                     ci.lb = study.results[[i]]$`95%lb`,
                                     ci.ub = study.results[[i]]$`95%ub`,
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
    for (i in 1:length(study.results)){
        true.C = study.results[[i]][4,2]
        if (true.C >= study.results[[i]][1,1] & true.C <= study.results[[i]][1,3]){
            cr.regular = cr.regular + 1
        }
        if (true.C >= study.results[[i]][2,1] & true.C <= study.results[[i]][2,3]){
            cr.logit = cr.logit + 1
        }
        if (true.C >= study.results[[i]][3,1] & true.C <= study.results[[i]][3,3]){
            cr.arcsine = cr.arcsine + 1
        }
    }
    # calculate proportion of times the c-statistic lies within the boostrap intervals
    cr.logit = cr.logit/length(study.results)
    cr.regular = cr.regular/length(study.results)
    cr.arcsine = cr.arcsine/length(study.results)
    return(list('regular.pooled'= cr.regular, 'logit.transformed' = cr.logit,
                'arcsine.transformed' = cr.arcsine))
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
                        "regular.pool", "pooled_se", "func", "arcsine.pool",
                        'n.imp','nr_boots'), 
                  envir=environment())
    # apply parallel computing
    study.results = parLapply(cl, sims, function(x, method=mech, on=var){
        # determine "true" c-statistic
        set.seed(11)
        true = boot::boot(x, func, R=50)
        #true.C.CI = boot::boot.ci(true, type = 'perc')$percent[c(4,5)]
        #true.C = c(true.C.CI[1], true$t0, true.C.CI[2])
        true.C = rep(true$t0, 3)
        set.seed(11)
        # create missing data according to user-specified mechanism
        missingdata = amputation(x, method, on)
        missingdata$V3 = as.factor(missingdata$V3) # change outcome to factor
        # Impute
        set.seed(11)
        invisible(capture.output(imp <- mice::mice(missingdata, n.imp, method = c('pmm', 'pmm', 'logreg'), seed = 5)))
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

