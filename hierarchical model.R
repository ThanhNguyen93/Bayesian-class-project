
library(LearnBayes)
data(hearttransplants)
attach(hearttransplants)
View(hearttransplants)

#PRIOR ~gamma(alpha, beta) has parameter alpha and beta
alpha = 16; beta = 15174
#prior mean
lam = alpha/beta
y = 0:10 #prior data from 10 hospitals

#SAMPLING from hospital A
yobs = 1; ex=66

#PREDICT mortality of 10 hospitals A
#predictive density = (sampling*prior)/posterior
py=dpois(y, lam*ex) *dgamma(lam, shape= alpha, rate = beta)/dgamma(lam, 
                                                                   shape = alpha + y, rate = beta + ex )
cbind(y, round(py, 3))

#The posterior density of λ can be summarized by simulating 1000 values 
#from the gamma density
lambdaA = rgamma(1000, shape = alpha + yobs, rate = beta+ex)

#consider another hopsital B with yobs=4, exposure e = 1767
ex = 1767; yobs=4
y=0:10
py=dpois(y, lam*ex)*dgamma(lam, shape = alpha, rate = beta)/dgamma(lam, 
                                                                   shape=alpha+y, rate = beta+ex)
cbind(y, round(py, 3))

#posterior mean
lambdaB= rgamma(1000, shape=alpha+yobs, rate =beta+ex)

#graph
lambda=seq(0, max(c(lambdaA, lambdaB)), length=500)
par(mfrow=c(2, 1))
hist(lambdaA, freq=FALSE, main=" ", ylim = c(0, 1500))
lines(lambda, dgamma(lambda, shape = alpha, rate = beta))
#hospital A with relatively little experience in surgeries, the prior information is significant
#and the posterior distribution resembles the prior distribution

hist(lambdaB, freq= FALSE, main=" ", ylim=c(0, 1500))
lines(lambda, dgamma(lambda, shape=alpha, rate = beta))
#for hospital B with many surgeries, the prior information is less influential 
#and the posterior distribution resembles the likelihood function



####chap 7: HIERARCHIAL MODEL
#y_i/e_i : # of deaths per unit exposure
#PLOT ratios {yi/ei} against the logarithms of the exposures {log(ei)} for all hospitals 
#where each point is labeled by the number of observed deaths yi.
plot(log(e), y/e, xlim=c(6,9.7), xlab="log(e)", ylab="y/e")
text(log(e),y/e,labels=as.character(y),pos=4)

sum(y)
sum(e)

#we simulate 1000 values from the posterior predictive density of y*_94.
#first, simulate 1000 draws of the posterior density of λ
lambda=rgamma(1000,shape=277,rate=294681)
# then simulate draws of y*_94 from a Poisson distribution with mean e_94*λ.
ys94=rpois(1000,e[94]*lambda)

# displays a histogram of this posterior predictive distribution
hist(ys94,breaks=seq(0.5,max(ys94)+0.5))
#actual number of transplant deaths y94 is shown by a vertical line
lines(c(y[94],y[94]),c(0,120),lwd=3)
#--> observed yj is in the tail portion of the distribution, it seems inconsistent with the fitted model – 
#it suggests that this hospital actually has a higher true mortality rate than estimated from this equal-rates model.


#next, check consistency of observed yi with its posterior predictive for ALL hospitals
lambda=rgamma(1000,shape=277,rate=294681)
prob.out=function(i)
{
  ysi=rpois(1000,e[i]*lambda)
  pleft=sum(ysi<=y[i])/1000
  pright=sum(ysi>=y[i])/1000
  min(pleft,pright)
}
pout=sapply(1:94,prob.out)
plot(log(e),pout,ylab="Prob(extreme)")
#--> Note that a number of these tail probabilities appear small (15 are smaller than 0.10) 
#which means that the “equal rates” model is inadequate for explaining the distribution of mortality rates 
#for the group of 94 hospitals

