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


## In Bayesian inference, we use probability directly to characterize
## uncertainty about parameters whereas in classical inference probability is used to characterize operating
## properties of statistical procedures.

## Bayesian inference is characterized by two features: 
## 1) posterior inference
## 2) regarding all parameters as realizations of random variables

## Bayes rule 
## [z|y] = [y|z][z]\[y]
## posterior = likelihood*prior/marginalization

## [z|y] = conditional probability distribution of z given y (the probability of z being true given y)
## [y|z] = conditional probability distribution of y given z (the probability of y being true given z)
## [z] = marginal distribution of z (the probability of z being true)
## [y] = marginal distribution of y (the probability of y being true)

## we can think of this as a way of updating our prior knowledge about an unknown
## quantity (expressed by [z]) with information obtained from the data ([y|z]).

## for example, in occupancy modeling terms:
## Pr(present|not detected) = Pr(not detected|present)*Pr(present) \ Pr(not detected)
## the probability a species is present given the absense of detection (the posterior)
## equals the probability that a species is not detected given being present (the likelihood)
## multiplied by the probability that a species is present (the prior)
## divided by the probability that a species is not detected (marginalization)

## another example in our modeling context:
## [model|data] = [data|model][model] \ [data]

## Regarding choosing priors:
## majority of time, chosen to express a lack of specific information to 'let the data speak'
## however can be chosen based on information from a comparable previous study, i.e. sequentially across years
## informative priors can be useful if you have a parameter that is technically undefinable
## according to Kery: we generally recommend the use of priors that are meant to reflect a lack of information

## TODO QUESTION: For our lab's modeling approach, what type of priors do we typically use?

## 2.6 Basic Markov Chain Monte Carlo (MCMC)

## a class of methods for drawing random samples (i.e., simulating) from the target posterior distribution (or any distribution for that matter).

## any particular MCMC algorithm needs to run for some number of iterations before the
## Markov chain “converges” to the target posterior distribution. Samples from this early, transitional
## period, called the “burn-in” or “warm-up,” are discarded and not used in characterizing posterior
## summaries.

## 2.6.2 ILLUSTRATION: USING MH FOR A BINOMIAL MODEL
# Simulate data
set.seed(2016)
y <- rbinom(2, size=10, p = 0.5)
# Define the joint distribution ([ likelihood) which we will maximize
jointdis <- function(data, J, p){
  prod(dbinom(data, size=J, p=p))
}
# Posterior is proportional to likelihood times prior
posterior <- function(p, data, J, a, b){
  prod(dbinom(data, size=J, p=p)) * dbeta(p, a, b)
}

# Do 100,000 MCMC iterations using Metropolis algorithm
# Assume vague prior which is beta(1,1) [ Unif(0,1)
mcmc.iters <- 100000
out <- rep(NA, mcmc.iters)

# Starting value
p <- 0.2
# Begin the MCMC loop
for(i in 1:mcmc.iters){
  # Use a uniform candidate generator (not efficient)
  p.cand <- runif(1, 0, 1)
  # Alternative: random walk proposal
  # p.cand <- rnorm(1, p, 0.05) # Need to reject if > 1 or < 0
  # if(p.cand < 0 | p.cand > 1 ) next
  r <- posterior(p=p.cand, y, J=10, a=1, b=1) / posterior(p=p, y, J=10, a=1, b=1)
  # Generate a uniform r.v. and compare with "r", this imposes the
  # correct probability of acceptance
  if(runif(1) < r)
    p <- p.cand
  # Save the current value of p
  out[i] <- p
}

mean(out)

sd(out)

quantile(out, c(0.025, 0.975))

# Evaluate likelihood for a grid of values of p
p.grid <- seq(0.1, 0.9, , 200)
likelihood <- rep(NA, 200)
for(i in 1:200){
  likelihood[i] <- jointdis(y, J=10, p=p.grid[i])
}

par(mfrow=c(2,1), mar = c(5,5,3,2))
plot(p.grid, likelihood, xlab="", ylab="Likelihood", xlim=c(0,1), ty = "l", main =
       "Likelihood function")
p.hat <- p.grid[likelihood == max(likelihood)]
abline(v = p.hat)
text(p.hat, 0.005, paste("MLE = ", round(p.hat, 3), sep= ""))
plot(density(out), xlim=c(0,1), main = "Posterior distribution", xlab = "p",
     ylab = "Posterior")
p.mean <- mean(out)
abline(v = p.mean)
text(p.mean, 0.5, paste("Post. mean = ", round(p.mean, 3),sep=" "))

## 2.6.3 METROPOLIS ALGORITHM FOR MULTIPARAMETER MODELS

log.posterior <- function(beta0, beta1, z, vegHt){
  # Note: "z" and "vegHt" must be input
  loglike <- -1 * negLogLike(c(beta0, beta1), z, vegHt)
  logprior <- dnorm(c(beta0, beta1), 0, 10, log=TRUE)
  return(loglike + logprior[1] + logprior[2])
}

niter <- 50000
out <- matrix(NA, niter, 2, dimnames = list(NULL, c("beta0", "beta1")))
# Initialize parameters
beta0 <- rnorm(1)
beta1 <- rnorm(1)
# Current value of the log(posterior)
logpost.curr <- log.posterior(beta0, beta1, z, vegHt)
# Run MCMC algorithm
for(i in 1:niter){
  if(i %% 1000 == 0) # report progress
    cat("iter", i, "yn")
  # Update intercept (beta0)
  # Propose candidate values of beta
  # If the proposal was not symmetric, would be Metrop-*Hastings*
  beta0.cand <- rnorm(1, beta0, 0.3) # 0.3 is tuning parameter
  # Evaluate the log(posterior)
  logpost.cand <- log.posterior(beta0.cand, beta1, z, vegHt)
  # Compute Metropolis acceptance probability, r
  r <- exp(logpost.cand - logpost.curr)
  # Keep candidate if it meets criterion (u < r)
  if(runif(1) < r){
    beta0 <- beta0.cand
    logpost.curr <- logpost.cand
  }
  # Update slope (beta1)
  beta1.cand <- rnorm(1, beta1, 0.3) # 0.3 is tuning parameter
  # Evaluate the log(posterior)
  logpost.cand <- log.posterior(beta0, beta1.cand, z, vegHt)
  # Compute Metropolis acceptance probability
  r <- exp(logpost.cand - logpost.curr)
  # Keep candidate if it meets criterion (u < r)
  if(runif(1) < r){
    beta1 <- beta1.cand
    logpost.curr <- logpost.cand
  }
  out[i,] <- c(beta0, beta1) # Save samples for iteration
}
# Plot
layout(rbind(c(1,1),
             c(2,2),
             c(3,4)), respect = T) # <- play with these settings
par(oma=c(0,0,0,0),mar=c(5,4,1,1))
plot(out[,1], type="l", xlab="Iteration", ylab="beta0")
plot(out[,2], type="l", xlab="Iteration", ylab="beta1")
plot(density(out[,1]), xlab="beta0", main="")
plot(density(out[,2]), xlab="beta1", main="")

## 2.6.7 BAYESIAN ANALYSIS OF HMs
psi <- plogis(-3 + 2*vegHt) # occupancy probability
z <- rbinom(M, 1, psi) # realised presence/absence
p <- 0.6 # detection probability
J <- 3 # sample each site 3 times
y <-rbinom(M, J, p*z) # observed detection frequency
# Number of MCMC iterations to to
niter <- 50000
# Matrix to hold the simulated values
out <- matrix(NA, niter, 3, dimnames = list(NULL, c("beta0", "beta1", "p")))
# Initialize parameters, likelihood, and priors
starting.values <- c(beta0=0, beta1=0)
beta0 <- starting.values[1]
beta1 <- starting.values[2]
z <-ifelse(y>0, 1, 0)
p <- 0.2

# NOTE: using logistic reg. likelihood function here (defined previously)
loglike <- -1*negLogLike(c(beta0, beta1), z, vegHt)
logprior <- dnorm(c(beta0, beta1), 0, 10, log=TRUE)
# Run MCMC algorithm
for(i in 1:niter) {
  if(i %% 1000 ==0 ) # report progress
    cat("iter", i, "yn")
  # PART 1 of algorithm -- same as before
  # Update intercept (beta0)
  # propose candidate values of beta
  beta0.cand <- rnorm(1, beta0, 0.3) # 0.3 is tuning parameter
  # evaluate likelihood and priors for candidates
  loglike.cand <- -1*negLogLike(c(beta0.cand, beta1), z, vegHt)
  logprior.cand <- dnorm(beta0.cand, 0, 10, log=TRUE)
  # Compute Metropolis acceptance probability (r)
  r <- exp((loglike.cand+logprior.cand) - (loglike + logprior[1]))
  # Keep candidate if it meets the criterion
  if(runif(1) < r){
    beta0 <- beta0.cand
    loglike <- loglike.cand
    logprior[1] <- logprior.cand
  }
  # Update slope (beta1)
  beta1.cand <- rnorm(1, beta1, 0.3) # 0.3 is tuning parameter
  # evaluate likelihood and priors for candidates
  loglike.cand <- -1*negLogLike(c(beta0,beta1.cand), z, vegHt)
  logprior.cand <- dnorm(beta1.cand, 0, 10, log=TRUE)
  # Compute Metropolis acceptance probability r
  r <- exp((loglike.cand+logprior.cand) - (loglike + logprior[2]))
  # Keep the candidates if they meet the criterion
  if(runif(1) < r) {
    beta1 <- beta1.cand
    loglike <- loglike.cand
    logprior[2] <- logprior.cand
  }
  # Part 2 of the algorithm
  # update z. Note we only need to update z if y=0.
  # The full conditional has known form
  psi <- plogis(beta0 + beta1 * vegHt)
  psi.cond <- dbinom(0,J,p) * psi /(dbinom(0, J, p) * psi + (1-psi))
  z[y==0] <- rbinom(sum(y==0), 1, psi.cond[y==0])
  loglike <- -1 * negLogLike(c(beta0, beta1), z, vegHt)
  
  # Part 3: update p
  ## The commented code will update p using Metropolis
  ## loglike.p <- sum(log(dbinom(y[z==1],J,p)))
  ## p.cand <- runif(1, 0, 1)
  ## loglike.p.cand <- sum(log(dbinom(y[z==1], J, p.cand)))
  ## if(runif(1) < exp(loglike.p.cand-loglike.p))
  ## p <- p.cand
  ## This bit draws p directly from its full conditional
  p <- rbeta(1, 1+ sum(y), sum(z)*J +1 - sum(y) )
  # Save MCMC samples
  out[i,] <- c(beta0,beta1,p)
}
# Plot bivariate representation of joint posterior
pairs(out)
# Trace/history plots for each parameter (Fig. 2.8)
op <- par(mfrow=c(2,1))
plot(out[,1], type="l", xlab="Iteration", ylab="beta0")
abline(h = mean(out[,1]), col = "blue", lwd = 2)
abline(h = -3, col = "red", lwd = 2)
plot(out[,2], type="l", xlab="Iteration", ylab="beta1")
abline(h = mean(out[,2]), col = "blue", lwd = 2)
abline(h = 2, col = "red", lwd = 2)


## 2.7 Model Selection and Averaging


## 2.8.1 PARAMETRIC BOOTSTRAPPING EXAMPLE
sim.data <- function(beta0 = -3, beta1 = 2, p = 0.6, x=NULL){
  # Function allows input of covariate "x", or simulates new
  M <- 100
  if(is.null(x))
    vegHt <- runif(M, 1, 3) # uniform from 1 to 3
  # Suppose that occupancy probability increases with vegHt
  # The relationship is described (default) by an intercept of -3 and
  # a slope parameter of 2 on the logit scale
  # plogis is the inverse-logit (constrains us back to the [0-1] scale)
  psi <- plogis(beta0 + beta1*vegHt)
  # Now we simulated true presence/absence for 100 sites
  z <- rbinom(M, 1, psi)
  # Now generate observations
  J <- 3 # sample each site 3 times
  y <- rbinom(M,J,p*z)
  list(y=y, J=J, vegHt=vegHt)
}

# This is the negative log-likelihood based on the marginal distribution
# of y. It is the pmf of a zero-inflated binomial random variable.
#
negLogLikeocc <- function(beta, y, x, J) {

  beta0 <- beta[1]
  beta1 <- beta[2]
  p <- plogis(beta[3])
  psi <- plogis(beta0 + beta1*x)
  marg.likelihood <- dbinom(y, J, p) * psi + ifelse(y==0, 1, 0) * (1-psi)
  return(-sum(log(marg.likelihood)))
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


data <- sim.data() # Generate a data set
# Let's minimize the negative log-likelihood
starting.values <- c(beta0=0, beta1=0, logitp=0)
opt.out <- optim(starting.values, negLogLikeocc, y=data$y, x=data$vegHt,J=data$J,
                 hessian=TRUE)
(mles <- opt.out$par)
# Make a table with estimates, SEs, and 95% CI
mle.table <- data.frame(Est=mles,
                        SE = sqrt(diag(solve(opt.out$hessian))))
mle.table$lower <- mle.table$Est - 1.96*mle.table$SE
mle.table$upper <- mle.table$Est + 1.96*mle.table$SE
mle.table

# Define a fit statistic
fitstat <- function(y, Ey){
  sum((sqrt(y) - sqrt(Ey)))
}
# Compute it for the observed data
T.obs <- fitstat(y, J*plogis(mles[1] + mles[2]*vegHt)*plogis(mles[3]))
# Get bootstrap distribution of fit statistic
T.boot <- rep(NA, 100)
for(i in 1:100){
  # Simulate a new data set and extract the elements. Note we use
  # the previously simulated "vegHt" covariate
  data <- sim.data(beta0=mles[1],beta1=mles[2],p=plogis(mles[3]),x=vegHt)
  # Next we fit the model
  starting.values <- c(0,0,0)
  opt.out <- optim(starting.values, negLogLikeocc, y=data$y, x= data$vegHt, J=data$J,
                   hessian=TRUE)
  (parms <- opt.out$par)
  # Obtain the fit statistic
  T.boot[i]<- fitstat(y, J*plogis(parms[1] + parms[2]*vegHt)*plogis(parms[3]) )

}

(T.obs)

summary(T.boot)



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