## Chapter 2 Kery
## WHAT ARE HIERARCHICAL MODELS AND HOW DO WE ANALYZE THEM?

## 2.1 Introduction

## 2.2 Random Variables, Probability Density Functions, Statistical Models, Probability, and Statistical Inference

## 2.3 Hierarchical Models (HMs)

## 2.4 Classical Inference Based on Likelihood

## 2.5 Bayesian Inference

## 2.6 Basic Markov Chain Monte Carlo (MCMC)

## 2.7 Model Selection and Averaging

## 2.8 Assessment of Model Fit

## 2.9 Summary and Outlook

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
beta0 beta1
-2.539793 1.617025
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