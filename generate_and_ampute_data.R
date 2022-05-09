# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")

# import packages
library(lattice)
library(mice)
library(tidyverse)
source('FUNCTIONS.R')

## GENERATE DATA ----
set.seed(11)
testdata = simulate.multivariate(nr_conditions = 'linear', threshold = 0.6, prob = 0.9)
# descriptives
summary(testdata)
#visualize 
plot(testdata$V1, testdata$V2, col = ifelse(testdata$V3==1,'green', 'red'))

# save as csv
## Save complete  dataframe as CSV file
write.csv(testdata,"C:\\Users\\surface\\Documents\\Applied Data Science\\Thesis\\complete_data.csv", row.names = FALSE)

# AMPUTA DATA ----
#                     V1 V2 V3   V1 V2 V3    V1 V2 V3    V1 V2 V3
pattern_list = list(c(0, 1, 1), c(1, 0, 1), c(1, 0, 1), c(0, 0, 0 ))
prop_list = c(0.8,0.8,0.8,0.8) # proportion missing values according to some missing mechanism
mech_list = c('MAR', 'MAR', 'MNAR', 'MCAR') # missing mechanisms underlying missingess
probabilities = c(0.4, 0.4, 0.1, 0.1) # relative occurrence of the missing mechanism in the data

# apply amputations function
set.seed(11)
ampdata = amputation(testdata, pattern_list, prop_list, mech_list, probabilities)
summary(ampdata)

#save missing data
write.csv(ampdata,"C:\\Users\\surface\\Documents\\Applied Data Science\\Thesis\\missingdata.csv", row.names = FALSE)

# Visualize ----
ggpubr::ggarrange(
    testdata %>% ggplot(aes(V1,V2, color = factor(V3))) + geom_point() + geom_smooth(method = 'lm'),
    na.omit(ampdata) %>% ggplot(aes(V1,V2, color = factor(V3))) + geom_point() + geom_smooth(method = 'lm')
)

