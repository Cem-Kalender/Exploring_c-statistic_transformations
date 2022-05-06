# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")

# import packages
library(lattice)
library(mice)
library(tidyverse)

## GENERATE DATA ----


simulate.multivariate <- function(n=10000, mu=c(22, 12, 0), sigma = c(1.0, 0.3, 0.5, 
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
    else {
        stop('Number of conditions not correctly specified')
    }
    return(data)
    }

set.seed(11)
testdata = simulate.multivariate()
# descriptives
summary(testdata)
#visualize 
plot(testdata$V1, testdata$V2, col = ifelse(testdata$V3==1,'green', 'red'))

# save as csv
## Save complete  dataframe as CSV file
write.csv(testdata,"C:\\Users\\surface\\Documents\\Applied Data Science\\Thesis\\complete_data.csv", row.names = FALSE)


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

# amputation function inputs
#                      V1 V2 V3   V1 V2 V3    V1 V2 V3    V1 V2 V3
pattern_list = list(c(0, 1, 1), c(1, 0, 1), c(1, 0, 1), c(0, 0, 0 ))
prop_list = c(0.8,0.8,0.8,0.8) # proportion missing values according to some missing mechanism
mech_list = c('MAR', 'MAR', 'MNAR', 'MCAR') # missing mechanisms underlying missingess
probabilities = c(0.4, 0.4, 0.1, 0.1) # relative occurrence of the missing mechanism in the data

# now test amputations function
set.seed(11)
ampdata = amputation(testdata, pattern_list, prop_list, mech_list, probabilities)
summary(ampdata)

# Visualize
ggpubr::ggarrange(
    testdata %>% ggplot(aes(V1,V2, color = factor(V3))) + geom_point() + geom_smooth(method = 'lm'),
    na.omit(ampdata) %>% ggplot(aes(V1,V2, color = factor(V3))) + geom_point() + geom_smooth(method = 'lm')
)


#save missing data
write.csv(ampdata,"C:\\Users\\surface\\Documents\\Applied Data Science\\Thesis\\missingdata.csv", row.names = FALSE)
