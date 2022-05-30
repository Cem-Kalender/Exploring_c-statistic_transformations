logitTransform <- function(p) { log(p/(1-p)) }
p <- seq(0.001, 0.999, 0.001)
pLogit <- logitTransform(p)
plot(p, pLogit, type='l', lwd=2, col='red', las=1, xlab='p', ylab='logit(p)')

asinTransform <- function(p) { asin(sqrt(p)) }
pAsin <- asinTransform(p)
plot(p, pAsin, type='l', lwd=2, col='blue', las=1, xlab='p', ylab='arcsine(p)')

rangeScale <- function(x) { (x-min(x)) / (max(x)-min(x)) }

pAsin.scaled <- rangeScale(pAsin)
pLogit.scaled <- rangeScale(pLogit)

plot(p, pAsin.scaled, las=1, type='l', lwd=2, col='blue', xlab='c-statistic', ylab='c-statistic transformed')
points(p, pLogit.scaled, type='l', lwd=2, col='red')
text(0.8, 0.8, 'asin', col='blue')
text(0.8, 0.5, 'logit', col='red')

