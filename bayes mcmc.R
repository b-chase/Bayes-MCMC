require(tidyverse)
library(coda)

n <- 500
x <- runif(n, -10, 30)
z <- runif(n, -30, -12)
b1 <- 27
b1.alt <- -5
b2 <- 0.04
a <- -2

prob.fail <- 0.0
y <- ifelse(runif(n)>prob.fail,b1*x,b1.alt*x) + b2*z + a + rnorm(n, 0, 1)

data <- tibble(y,x,z)

lm(y~x+z, data)

metro_mc <- function(data, iter = 10000, burn = 500, starting=NULL, proposals="rwalk", priors="vague") {
  # data is a table with dependent variable in first column
  # mmc_spec is a list of two functions, priors and LH
  # prior is a function for drawing potential values of the parameters
  # LH is the likehood, a function that gives p density for a given error
  # iter is number of MC steps to take
  # burn is burn-in, number of initial steps to ignore
  
  k <- ncol(data) #n-1 independent variables plus a constant
  var.names <- c(names(data[,2:k]), "intercept", "nv")
  # k+1 terms if we count nv from t-distributed error (can also use precision or somethign else)
  
  if (is.null(starting)) {starting <- c(rep(0,k), 5)}
  names(starting) <- var.names
  if (is.character(proposals)) {
    proposals <- function(old.props) {new <- c(old.props[1:k] + rnorm(k,-10,10), rgamma(1, 1, 10))}
  }
  if (is.character(priors)) {
    priors <- function(params) {c(rep(1,k), 1/params[k+1])}
  }
  
  data <- mutate(data, const.=rep(1,nrow(data))) #necessary for next step
  
  likelihood <- function(inputs, params) {
    k <- length(params) - 1
    nv <- params[k+1]
    deviates <- apply(inputs, 1, function(pnt) 
      {pnt[1] - params[1:k] %*% pnt[2:(k+1)]})
    #print(summary(deviates))
    #print(nv)
    LLH <- sum(log(dt(deviates, nv)))
    LLH
  }
  
  #list of all parameters that show up in chain
  param.list <- as.matrix(rbind(starting)) 
  old.LLH <- likelihood(data, starting)
  #let's chain
  accepts <- 0
  
  for (i in seq(iter)) {
    # calculate prior densities
    old.params <- param.list[i,]
    #print(old.params)
    proposed.params <- proposals(old.params)
    #names(proposed.params) <- var.names
    #print(proposed.params)
    old.prior.prob <- prod(priors(old.params))
    new.prior.prob <- prod(priors(proposed.params))
    new.LLH <- likelihood(data, proposed.params)
    alpha.numer <- new.LLH + log(new.prior.prob)
    alpha.denom <- old.LLH + log(old.prior.prob)
    #print(alpha.numer)
    #print(alpha.denom)
    alpha <- alpha.numer - alpha.denom
    if (is.na(alpha)) {
      # output error information
      print(proposed.params)
      print(alpha.numer)
      print(old.params)
      print(alpha.denom)
    }
    #print(alpha)
    new.row <- c()
    if (log(runif(1)) < alpha) {
      new.row <- proposed.params
      accepts <- accepts + 1
      old.LLH <- new.LLH
    } else {
      new.row <- old.params
    }
    #print(new.row)
    param.list <- rbind(param.list, new.row)
    #print("Full list:")
    #print(param.list)
    if (i %% 50 == 0) {
      cat(paste0("\r", "Number of Accepts / Proposals: ", accepts, " / ", i))
    }
  }
  cat("\nFinished\n")
  as_tibble(param.list[(burn+1):(iter+1),])
}


auto_mc <- function(data, batch.sizes=5000) {
  k <- ncol(data) #number of variables, not including intercept
  print("First running a test chain of 2000 samples (burning first 500):")
  test.chain <- metro_mc(data, iter=2000, burn=500)
  summary(test.chain)
  test.accepts <- nrow(unique(test.chain[,1]))
  test.accepts
  test.means <- apply(test.chain, 2, mean)
  test.var <- sapply(test.means, function(t) max(abs(t)/1.6, 1))
  if (test.accepts >= 3) {
    test.var <- apply(test.chain, 2, function(m) {10*(max(m)-min(m))/sqrt(3)})
    # i know this is messy, but using ifelse wasn't working, will try fixing it later
  }
  print(test.means)
  print(test.var)
  
  zone_in <- function(data.tib, last.chain, means, varis, batches=5, batch.size=5000, batch.burn=500) {
    k <- ncol(data.tib)
    if (batches == 0) {
      return(last.chain)
    } else {
      cat(paste0("\nBatches remaining: ", batches, "\n"))
      chain.proposals <- function(old.params) {c(rnorm(k,means[1:k],2*varis[1:k]), 
                                                 rgamma(1,means[k+1]^2/varis[k+1],varis[k+1]/means[k+1]))}
      
      new.chain <- metro_mc(data, iter=batch.size, burn=batch.burn, proposals=chain.proposals, starting=means)
      
      chain.accepts <- nrow(unique(new.chain[,1]))
      chain.means <- apply(new.chain, 2, mean)
      chain.varis <- sapply(chain.means, function(t) max(abs(t)/2, 1))
      if (chain.accepts >= 3) {
        chain.varis <- apply(new.chain, 2, function(m) {5*(max(m)-min(m))/sqrt(3)})
        # i know this is messy, but using ifelse() wasn't working, will try fixing it later
      }
      print(rbind(chain.means, chain.varis))
      return(do.call(zone_in, args=list(data.tib, new.chain, chain.means, chain.varis, batches-1, batch.size, batch.burn)))
    }
  }
  
  last.chain <- zone_in(data, test.chain, test.means, test.var, batch.size = batch.sizes, batch.burn = 0)
  par(mfcol=c(2,ceiling(k/2)+1))
  traceplot(mcmc(last.chain))  
  
  last.chain
}

my.chain <- auto_mc(data)
summary(my.chain)
