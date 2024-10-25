## Chapter 2 Kery
## WHAT ARE HIERARCHICAL MODELS AND HOW DO WE ANALYZE THEM?

## 2.1 Introduction

## 2.2 Random Variables, Probability Density Functions, Statistical Models, Probability, and Statistical Inference

## 2.3 Hierarchical Models (HMs)

## 2.4 Classical Inference Based on Likelihood

## 2.5 Bayesian Inference

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