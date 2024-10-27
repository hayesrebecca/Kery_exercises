## Chapter 2 Kery
## WHAT ARE HIERARCHICAL MODELS AND HOW DO WE ANALYZE THEM?

## 2.1 Introduction

## Hierarchical model <- set of models that are conditionally related.

## Probability distribution of one random variable depends on the other random 
## variable. Models contain one component for the observations and one or more 
## components that describe the variables or outcomes of ecological process.

## 2.2 Random Variables, Probability Density Functions, Statistical Models, Probability, and Statistical Inference

## random variables <- variables for which possible values are governed by
## probability distributions.

## Two types of probability distributions:
## probability density function(pdf) for continuous variables
## probability mass function (pmf) for discrete variables
## They depend on one or more quatities called parameters

## Example of a PMf (counts of observed peregrine falcons)
## y = {0,1,2,3,4,5}, N = 5, p = 0.2
## Probability that we detect a falcon 0-5 times
dbinom(0:5, size = 5, prob = 0.2)

## Example of a PDF (height of adult males)
## mean = 190 cm, standard dev = 10 cm
##Proportion of the population shorter than 200cm
pnorm(200, mean = 190, sd = 10)
# Expected proportion of males with height between 180-200
pnorm(200, mean = 190, sd = 10)- pnorm(180, mean = 190, sd = 10)

## Statistical model <- observed response is treated as a random variable with
## an assume probability distribution (pdf) and this response can be predicted in
## average but there is variability that is often caused by other variables(covariates). 

## Joint probability distribution

## Probability of simultaneously observing all unique combinations of Y and X. 
## Written as [Y = y, X = x]
## Example: 
## Sets up the observations
Y <- 0:5 #Possible values of Y (# surveys with peregrine sightings)
X <- 0:5 #Possible values of X (# fledged young)
p <- plogis(-1.2 + 2*X) # p as function of X
round(p, 2)
# Joint distribution [Y, X]
lambda <- 0.4
joint <- matrix(NA, length(Y), length(X))
rownames(joint) <- paste("y=", Y, sep="")
colnames(joint) <- paste("x=", X, sep="")
for(i in 1:length(Y)) {
  joint[,i] <- dbinom(Y, 5, p[i]) * dpois(X[i], lambda)
}
round(joint, 3)

## Marginal probability distribution

## For a variable Y its distribution averaged over all possible values of X. 
## If X is discrete we do the averaging, if X is continuous we integrate over 
## the range of X
## Example:
margX <- colSums(joint)
round(margX, 4)

margY <- rowSums(joint)
round(margY, 4)

## Conditional probability distribution
## Explicit description of the relationship between one random variables and another.
## Example:
YgivenX <- joint / matrix(margX, nrow(joint), ncol(joint), byrow=TRUE)
round(YgivenX, 2)

## 2.3 Hierarchical Models (HMs)

## Occupancy model for species distributions
## See page 11
## N-mixture models- used for animal abundance
## See page 12

## 2.4 Classical Inference Based on Likelihood

## "Statistical inference for any system is based on the conceptual view that data we observe (or may
## potentially observe) are outcomes of random variables having a distribution with parameters that we
## would like to learn about from data."

## Example:
# Simulate a covariate called vegHt for 100 sites
set.seed(2014) # Set seed so we all get the same values of vegHt
M <- 100 # Number of sites surveyed
vegHt <- runif(M, 1, 3) # uniform from 1 to 3

# Suppose that occupancy probability increases with vegHt
# The relationship is described by an intercept of -3 and
# a slope parameter of 2 on the logit scale
beta0 <- -3
beta1 <- 2
psi <- plogis(beta0 + beta1*vegHt) # apply inverse logit
# Now we go to 100 sites and observe presence or absence
z <- rbinom(M, 1, psi)

# Definition of negative log-likelihood.
negLogLike <- function(beta, y, x) {
  beta0 <- beta[1]
  beta1 <- beta[2]
  psi <- plogis(beta0 + beta1*x)
  likelihood <- psi^y * (1-psi)^(1-y) # same as next line:
  # likelihood <- dbinom(y, 1, psi)
  return(-sum(log(likelihood)))
}
# Look at (negative) log-likelihood for 2 parameter sets
negLogLike(c(0,0), y=z, x=vegHt)
negLogLike(c(-3,2), y=z, x=vegHt) # Lower is better!
# Let's minimize it formally by function minimisation
starting.values <- c(beta0=0, beta1=0)
opt.out <- optim(starting.values, negLogLike, y=z, x=vegHt, hessian=TRUE)
(mles <- opt.out$par) # MLEs are pretty close to truth

# Alternative 1: Brute-force grid search for MLEs
mat <- as.matrix(expand.grid(seq(-10,10,0.1), seq(-10,10,0.1)))
# above: Can vary resolution (e.g., from 0.1 to 0.01)
nll <- array(NA, dim = nrow(mat))
for (i in 1:nrow(mat)){
  nll[i] <- negLogLike(mat[i,], y = z, x = vegHt)
}
which(nll == min(nll))
mat[which(nll == min(nll)),]
# Produce a likelihood surface, shown in Fig. 2-2.
library(raster)
r <- rasterFromXYZ(data.frame(x = mat[,1], y = mat[,2], z = nll))
mapPalette <- colorRampPalette(rev(c("grey", "yellow", "red")))

plot(r, col = mapPalette(100), main = "Negative log-likelihood",
     xlab = "Intercept (beta0)", ylab = "Slope (beta1)")
contour(r, add = TRUE, levels = seq(50, 2000, 100))
# Alternative 2: Use canned R function glm as a shortcut
(fm <- glm(z ~ vegHt, family = binomial)$coef)
# Add 3 sets of MLEs into plot
# 1. Add MLE from function minimisation
points(mles[1], mles[2], pch = 1, lwd = 2)
abline(mles[2],0) # Put a line through the Slope value
lines(c(mles[1],mles[1]),c(-10,10))
# 2. Add MLE from grid search
points(mat[which(nll == min(nll)),1], mat[which(nll == min(nll)),2],
       pch = 1, lwd = 2)
# 3. Add MLE from glm function
points(fm[1], fm[2], pch = 1, lwd = 2)

Vc <- solve(opt.out$hessian) # Get variance-covariance matrix
ASE <- sqrt(diag(Vc)) # Extract asymptotic SEs
print(ASE)

# Compare to SEs reported by glm() function (output thinned)
summary(glm(z ~ vegHt, family = binomial))


## Bootstrapping
nboot<-1000 #Obtain1000bootstrapsamples
boot.out<-matrix(NA,nrow=nboot,ncol=3)
dimnames(boot.out)<-list(NULL,c("beta0","beta1","psi.bar"))
for(i in 1:1000){
  #Simulatedata
  psi<-plogis(mles[1]+mles[2]*vegHt)
  z<-rbinom(M,1,psi)
  #Fitmodel
  tmp<-optim(mles,negLogLike,y=z,x=vegHt,hessian=TRUE)$par
  psi.mean<-plogis(tmp[1]+tmp[2]*mean(vegHt))
  boot.out[i,]<-c(tmp,psi.mean)
}
SE.boot<-sqrt(apply(boot.out,2,var)) #GetbootstrapSE
names(SE.boot)<-c("beta0","beta1","psi.bar")
apply(boot.out,2,quantile,c(0.025,0.975))

SE.boot

## Discrete Random effects
set.seed(2014)
M <- 100
vegHt <- runif(M, 1, 3)
# number of sites
# uniform from 1 to 3
psi <- plogis(beta0 + beta1 * vegHt) # occupancy probability
z <- rbinom(M, 1, psi)
p <- 0.6
J <- 3
y <-rbinom(M, J, p*z)
# Define negative log-likelihood.
# realised presence/absence
# detection probability
# sample each site 3 times
# observed detection frequency
negLogLikeocc <- function(beta, y, x, J) {
  beta0 <- beta[1]
  beta1 <- beta[2]
  p <- plogis(beta[3])
  psi <- plogis(beta0 + beta1*x)
  marg.likelihood <- dbinom(y, J, p)*psi + ifelse(y==0,1,0)*(1-psi)
  return(-sum(log(marg.likelihood)))
}
starting.values <- c(beta0=0, beta1=0,logitp=0)
(opt.out <- optim(starting.values, negLogLikeocc, y=y, x=vegHt,J=J,
                  hessian=TRUE))
## To optain standard errors
sqrt(diag(solve(opt.out$hessian)))

## Continuous Latent Variable (random effect)
marg <- rep(NA, J+1)
for(j in 0:J){
  marg[j+1] <- integrate(
    function(x){
      dbinom(j, J, plogis(x)) * dnorm(x, mu, sigma)},
    lower=-Inf,upper=Inf)$value
  }

#nx[encounter frequencies, number inds. encountered 1, 2, ..., 14 times
nx <- c(34, 16, 10, 4, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0)
nind <- sum(nx)
# Number of individuals observed
J <- 14
# Number of sample occasions
# Model Mh likelihood
Mhlik <- function(parms){
  mu <- parms[1]
  sigma <- exp(parms[2])
  # n0 = number of UNobserved individuals: N = nind + n0
  n0 <- exp(parms[3])
  # Compute the marginal probabilities for each possible value j[0,..,14
  marg <- rep(NA,J+1)
  }
  
# The negative log likelihood involves combinatorial terms computed
# using lgamma()-1*(lgamma(n0 + nind + 1)- lgamma(n0 + 1) + sum(c(n0, nx) * log(marg)))

(tmp <- nlm(Mhlik, c(-1, 0, log(10)), hessian=TRUE))

(SE <- sqrt( (exp(tmp$estimate[3])^2)* diag(solve(tmp$hessian))[3] ) )

## 2.5 Bayesian Inference

## 2.6 Basic Markov Chain Monte Carlo (MCMC)

## 2.7 Model Selection and Averaging

## 2.8 Assessment of Model Fit

## 2.9 Summary and Outlook


## Exercises

# 1. You should be able to apply Bayes’ rule to the peregrine falcon example earlier in this chapter to
# compute the distribution of XrY. That is: how many fledged young are there, given that we have
# detected birds on 0, 1, 2, . , J visits?



# 2. For the bootstrap GoF analysis done on the occupancy model, we found that the model appears to
# fit the data well. Try fitting the wrong model; i.e., without the vegHt covariate, and see if the
# model fails the GoF test.



# 3. Where we talked about prior distributions and sensitivity we said “if we have a flat prior on
# logit(p) for some probability parameter p, this is very different from having a Uniform(0,1) prior
# on p.” Evaluate this by simulating data for different priors on logit(p) and then back-transforming
# the simulated values to see what the implied prior is for p. Find a normal prior for logit(p) that is
# approximately uniform on [0,1] for p.



# 4. Use a Metropolis or MH algorithm to simulate Normal(0,1) random variables if all you have
# access to is a uniform random number generator. Of course you also know the mathematical form
# for the normal pdf but you do not have a normal random number generator at hand. Verify that
# your simulated data have the required normal distribution by making a histogram and computing
# summary statistics.