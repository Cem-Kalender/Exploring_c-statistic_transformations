# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")

# import packages
library(lattice)
library(mice)
library(tidyverse)

## GENERATE DATA ----

set.seed(12)
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
testdata$V3 = ifelse(condition,rbinom(nrow(testdata[condition,]), size = 1, prob = 0.33),
                     ifelse(condition2,rbinom(nrow(testdata[condition2,]), size = 1, prob = 0.05), 0))
plot(testdata$V1, testdata$V2, col = ifelse(testdata$V3==1,'red', 'green'))

# descriptives
summary(testdata)
# save as csv
## Save complete  dataframe as CSV file
write.csv(testdata,"C:\\Users\\surface\\Documents\\Applied Data Science\\Thesis\\complete_data.csv", row.names = FALSE)


## AMPUTE WITH MICE ----

# ampute the complete data once for every mechanism
ampdata1 <- ampute(testdata, patterns = c(0, 1, 1), prop = 0.7, mech = "MAR")$amp
ampdata2 <- ampute(testdata, patterns = c(1, 0, 1), prop = 0.7, mech = "MNAR")$amp
ampdata3 <- ampute(testdata, patterns = c(0, 1, 1), prop = 0.7, mech = "MNAR")$amp
ampdata4 <- ampute(testdata, patterns = c(0, 0, 0), prop = 0.7, mech = "MCAR")$amp

# create a random allocation vector
# use the prob argument to specify how much of each mechanism should be created
# here, 0.5 of the missingness should be MAR and 0.5 should be MCAR
indices <- sample(x = c(1, 2, 3, 4), size = nrow(testdata), 
                  replace = TRUE, prob = c(0.6, 0.2, 0.1, 0.1))

# create an empty data matrix
# fill this matrix with values from either of the two amputed datasets
ampdata <- matrix(NA, nrow = nrow(testdata), ncol = ncol(testdata))
ampdata[indices == 1, ] <- as.matrix(ampdata1[indices == 1, ])
ampdata[indices == 2, ] <- as.matrix(ampdata2[indices == 2, ]) 
ampdata[indices == 3, ] <- as.matrix(ampdata3[indices == 3, ])
ampdata[indices == 4, ] <- as.matrix(ampdata4[indices == 4, ]) 

# store as df
ampdata = data.frame(ampdata)
# change colnames
colnames(ampdata) = c('V1', 'V2', 'V3')


# Visualize
ggpubr::ggarrange(
    testdata %>% ggplot(aes(V1,V2, color = factor(V3))) + geom_point() + geom_smooth(method = 'lm'),
    na.omit(ampdata) %>% ggplot(aes(V1,V2, color = factor(V3))) + geom_point() + geom_smooth(method = 'lm')
)


#save missing data
write.csv(ampdata,"C:\\Users\\surface\\Documents\\Applied Data Science\\Thesis\\missingdata.csv", row.names = FALSE)
