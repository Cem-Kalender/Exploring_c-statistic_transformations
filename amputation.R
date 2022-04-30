# Setting working directory and loading libraries----
setwd("~/Applied Data Science/Thesis")

# import packages
library(lattice)
library(mice)
library(tidyverse)

## simulate dataset----

set.seed(2022)
#options(scipen=999)

# randomly generate a dataset to be our complete dataset


testdata <- as.data.frame(MASS::mvrnorm(n = 10000, 
                                        mu = c(22, 12, 0),      # V1   V2   V3
                                        Sigma = matrix(data = c(1.0, 0.6, 0.4, 
                                                                0.6, 1.0, 0.5, 
                                                                0.4, 0.5, 1.0), 
                                                       nrow = 3, 
                                                       byrow = T)))
summary(testdata)

# People falling in a certain region have a a higher probabilityof Y = 1
condition = testdata$V1>median(testdata$V1) & testdata$V2>median(testdata$V2)
condition2 = testdata$V1<median(testdata$V1) & testdata$V2>median(testdata$V2)
# plot
plot(testdata$V1, testdata$V2, col = ifelse(condition, 'red', 'green'))

# generate binomial variable for V3
testdata$V3 = ifelse(condition,rbinom(nrow(testdata[condition,]), size = 1, prob = 0.5),
                     ifelse(condition2,rbinom(nrow(testdata[condition2,]), size = 1, prob = 0.1), 0))
plot(testdata$V1, testdata$V2, col = ifelse(testdata$V3==1,'red', 'green'))


# visualize missing mechanism plotting data points----
par(mfrow=c(1,4))
# no missingness
plot(testdata$V1, testdata$V2, main = 'Normal')
abline(lm(testdata$V2~testdata$V1), col = 'blue')
# MAR
plot(testdata$V1, testdata$V2, col = ifelse(testdata$V1>median(testdata$V1),'red','black'), main = 'Missingness: MAR')
abline(lm(testdata$V2[testdata$V1<median(testdata$V1)]~testdata$V1[testdata$V1<median(testdata$V1)], data = testdata), col = "blue")
# MNAR
plot(testdata$V1, testdata$V2, col = ifelse(testdata$V2>median(testdata$V2),'red','black'), main = 'Missingness: MNAR')
abline(lm(testdata$V2[testdata$V2<median(testdata$V2)]~testdata$V1[testdata$V2<median(testdata$V2)], data = testdata), col = "blue")

# create MAR+MNAR missingness mechanism and ampute ----

# first split the complete data, before amputing
set.seed(2022)
samp = sample(seq_len(floor(nrow(testdata)*0.8)))
test.data = testdata[-samp,] # this the test data that has not been amputed

# we will now ampute
missingdata = testdata[samp,]

# defining missing mechanisms
missing.mechanism1 = missingdata$V2<median(missingdata$V2) & missingdata$V1<median(missingdata$V1)
missing.mechanism2 = missingdata$V2>median(missingdata$V2) & missingdata$V1>median(missingdata$V1)
missinglist = ifelse(missing.mechanism1, rbinom(nrow(missingdata[missing.mechanism1,]), 1, 0.5),
                     ifelse(missing.mechanism2, rbinom(nrow(missingdata[missing.mechanism1,]), 1, 0.7),0))

# visualize
plot(missingdata$V1, missingdata$V2, col = ifelse(missinglist==1, 'red', 'green'),
     main = 'Incorporated missing mechanism')
abline(lm(testdata$V2~testdata$V1), col = "blue")

# Ampute (univariate missingness)
missingdata$V1[missinglist==1] = NA
abline(lm(missingdata$V2~missingdata$V1), col = "red")

sum(is.na(missingdata$V1)) # 3422 missing values created

# final visual comparison of complete datasets and amputed dataset
ggpubr::ggarrange(ggplot(testdata, aes(V1, V2, col = factor(V3))) + geom_point() + geom_smooth(col='black', method = 'lm'), 
                  ggplot(missingdata, aes(V1, V2, col = factor(V3))) + geom_point() + geom_smooth(col='black', method = 'lm'))

## Save amputed dataframe as CSV file
write.csv(missingdata,"C:\\Users\\surface\\Documents\\Applied Data Science\\Thesis\\missingdata.csv", row.names = FALSE)
