# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")

# import packages
library(lattice)
library(mice)
library(tidyverse)

## GENERATE DATA ----

set.seed(1202)
#options(scipen=999)

# randomly generate a dataset to be our complete dataset
testdata <- as.data.frame(MASS::mvrnorm(n = 10000, 
                                        mu = c(22, 12, 0),      # V1   V2   V3
                                        Sigma = matrix(data = c(1.0, 0.3, 0.5, 
                                                                0.3, 1.0, 0.5, 
                                                                0.5, 0.5, 1.0), 
                                                       nrow = 3, 
                                                       byrow = T)))

# People falling in a certain region have a a higher probabilityof Y = 1
condition = testdata$V1>median(testdata$V1) & testdata$V2>median(testdata$V2)
condition2 = testdata$V1<median(testdata$V1) & testdata$V2>median(testdata$V2)
# plot
plot(testdata$V1, testdata$V2, col = ifelse(condition, 'red', 'green'))

# generate binomial variable for V3
testdata$V3 = ifelse(condition,rbinom(nrow(testdata[condition,]), size = 1, prob = 0.5),
                     ifelse(condition2,rbinom(nrow(testdata[condition2,]), size = 1, prob = 0.1), 0))
plot(testdata$V1, testdata$V2, col = ifelse(testdata$V3==1,'red', 'green'))

# descriptives
summary(testdata)
# save as csv
## Save complete  dataframe as CSV file
write.csv(testdata,"C:\\Users\\surface\\Documents\\Applied Data Science\\Thesis\\complete_data.csv", row.names = FALSE)


## AMPUTE WITH MICE ----
result <- ampute(data = testdata)

is.na(missingdata) %>% sum()
md.pattern(result$amp)
# inspect missing patterns and insert another missingness pattern
mypatterns = result$patterns
mypatterns = rbind(mypatterns, c(0, 0, 1))
result <- ampute(testdata, patterns = mypatterns)

# the weights matrix can also be used to switch from MAR to MNAR missingness,
myweights <- result$weights
myweights[1, ] <- c(0, 0.8, 0.4)
myweights[3, ] <- c(3.0, 1.0, 0)
myweights
# ampute again
result <- ampute(testdata, weights = myweights, 
                 patterns = mypatterns, mech = "MNAR")

missingdata = result$amp

#save missing data
write.csv(missingdata,"C:\\Users\\surface\\Documents\\Applied Data Science\\Thesis\\univariate_missingdata.csv", row.names = FALSE)
