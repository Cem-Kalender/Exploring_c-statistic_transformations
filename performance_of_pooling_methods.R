# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")
# import packages
source('functions_simulation_study.R')
options(scipen=0)


# Read simulation study results from files
simulation_study = readRDS("simulation_study.RData")

##################################### COVERAGE ##############################
covered = coverage.process(simulation_study)
cov = coverage.rate(covered)
cov

# What are the means of theta hat?
t(sapply(c(2,1,3), function(i){
    
    theta_hat = sapply(simulation_study, function(x){
        sapply(x, function(z){
            z[i,1]})
        
    })
    theta_hat %>% colMeans()
}))


# how wide is the confidence interval per pooling method?
CI.broadness=lapply(covered, function(c){
    lapply(1:1000, function(x){
        coverage = c[[x]]
        logit.coverage = coverage[1,3]-coverage[1,1]
        regular.coverage = coverage[2,3]-coverage[2,1]
        arcsin.coverage = coverage[3,3]-coverage[3,1]
        rbind(regular.coverage, logit.coverage, arcsin.coverage)
    })
})

CI.width = sapply(CI.broadness, function(x){
    x=df.maker(x) # regular, logit, arcsine
    colMeans(x)
}) %>% data.frame()
colnames(CI.width) = c('0.7', '0.8', '0.9')
rownames(CI.width) = c('regular', 'logit', 'arcsine')
CI.width

# monte carlo error of coverage
monte.carlo.coverage = function(cover){
    return(sqrt((cover*(1-cover)/1000)))
}

coverage_df = cov %>% data.frame()
colnames(coverage_df) = c('true.c.0.7', 'true.c.0.8', 'true.c.0.9')
coverage_df
cbind(coverage_df[,1], monte.carlo.coverage(coverage_df[,1]),
      coverage_df[,2], monte.carlo.coverage(coverage_df[,2]),
      coverage_df[,3], monte.carlo.coverage(coverage_df[,3]))*100

################################# DISTRIBUTION OF POOLED ESTIMATES ###########################
# histograms of c-index estimates per value per sub-study
substudies = c('theta = 0.7', 'theta = 0.8', 'theta = 0.9')
method = c('logit', 'regular', 'arcsine')
C = c(0.7,0.8,0.9)

par(mfrow = c(1, 3))
sapply(1:3, function(g){
    sapply(c(2,1,3), function(h){
        c_statistic = sapply(1:1000, function(i){simulation_study[[g]][[i]][,1][h]})
        hist(c_statistic , breaks=20, main = '', cex.main = 1.5, xlab = '', ylab='', ylim=c(0,170))
        abline(v=mean(c_statistic ), lwd=2, col='yellow')
        abline(v=C[g], lwd=2, col='red')
    })
})

# determine skweness and kurtosis
par(mfrow = c(1, 3))
lapply(1:3, function(g){
    sapply(c(2,1,3), function(h){
        c_statistic = sapply(1:1000, function(i){simulation_study[[g]][[i]][,1][h]})
        c('skewness'=moments::skewness(c_statistic), 'kurtosis'=moments::kurtosis(c_statistic), 'sd'=sd(c_statistic))
    })
    
})

# distribution of pooled SE's ----

#Make scatterplots and determine Pearson's correlation between theta and se
par(mfrow = c(1, 3))
sapply(1:3, function(g){
    sapply(c(2,1,3), function(h){
        c_statistic = sapply(1:1000, function(i){simulation_study[[g]][[i]][,1][h]})
        SE = sapply(1:1000, function(i){simulation_study[[g]][[i]][,2][h]})
        plot(c_statistic,SE, xlab='Pooled estimates', ylab='Pooled SE')
        mtext(paste("Correlation:", round(cor(c_statistic, SE), 2)), cex=1, col='red')
        abline(lm(SE~c_statistic), lwd=3, col='red')
    })
    
})

# sub-study and pooling method-stratified regression model ----
se.vs.c = lapply(1:3, function(g){
    lapply(c(2,1,3), function(h){
        c_statistic = sapply(1:1000, function(i){simulation_study[[g]][[i]][,1][h]})
        SE = sapply(1:1000, function(i){simulation_study[[g]][[i]][,2][h]})
        lm(scale(SE)~scale(c_statistic)) %>% summary()
    })
    
})

# get coefficients and insert into 'mat'
count = 1
mat = NULL
for (i in 1:3){
    for(j in 1:3){
        mat[[count]]=se.vs.c[[i]][[j]]$coefficients['scale(c_statistic)',c('Estimate', 'Std. Error', 't value','Pr(>|t|)')]
        count = count + 1
    }
}
# covert 'mat' to dataframe
mat = data.frame(df.maker(mat))
colnames(mat) = c('Estimate', 'Std. Error', 't value','Pr(>|t|)')
rownames(mat) = paste0(rep(c('Regular', 'Logit', 'Arcsine'), 3),'/', rep(c('theta = 0.7', 'theta = 0.8', 'theta = 0.9'), each=3))
stargazer::stargazer(mat, type='text', summary=F)

######################################### BIAS #################################
biass = sapply(1:3, function(x){bias(simulation_study[[x]])}) %>% data.frame()
rownames(biass) = c('regular', 'logit', 'arcsine')
colnames(biass) = c('theta=0.7', 'theta=0.8', 'theta=0.9')

# monte carlo error of bias
monte.carlo.bias = function(theta){
    mc.se=sqrt((1/(1000*999))*sum((theta - mean(theta)**2)))
    return(mc.se)
}
# Get 'long' version of study results
simulation.long = lapply(c(2,1,3), function(z){
    sapply(1:3, function(x){
        lapply(1:1000, function(i){
            simulation_study[[x]][[i]][z]
        })
    })
    
})
# calculate monte carlo error of bias estimates
mcerrors=lapply(simulation.long, function(estimates){
    estimates = as.data.frame(estimates)
    apply(estimates, 2, function(x){monte.carlo.bias(theta=unlist(x))})
})
mcerrors=df.maker(mcerrors)
colnames(mcerrors) = c('theta=0.7', 'theta=0.8', 'theta=0.9')
rownames(mcerrors) = c('regular', 'logit', 'arcsine')

# inspect
biass
mcerrors

# Difference between theta hat and theta, per sub-study per pooling method
par(mfrow=c(1,1))
lapply(simulation_study, function(x){
    difference = lapply(x, function(i){
        
        theta = i[4,2]
        theta_hat_regular = i[2,1]
        theta_hat_logit = i[1,1]
        theta_hat_arcsine = i[3,1]
        
        r = theta_hat_regular - theta
        l = theta_hat_logit - theta
        a = theta_hat_arcsine - theta
        c(r,l,a)
    })
    sub.study = df.maker(difference)
    sapply(1:3, function(method){
        methods = c('Regular', 'Logit', 'Arcsine')
        # Create scatterplot
        plot(sub.study[,method], 1:1000, xlab='', cex.axis= 1.5, cex.lab = 2, ylab=methods[method], xlim=c(-0.1, 0.1))
        abline(v=mean(sub.study[,method]), lwd=3, col='red')
    })
})
