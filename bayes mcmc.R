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
  
  if (is.null(starting)) {starting <- c(rep(0,k), 5)}
  names(starting) <- var.names
  if (proposals[1] == "rwalk") {
    proposals <- sapply(seq(k), function(t) function(old) {old+rnorm(1,0,10)})
    proposals <- append(proposals, function(old) {runif(1,1,150)}) # last one is for the precision, tau
  }
  names(proposals) <- var.names
  if (priors[1] == "vague") {
    priors <- sapply(seq(k), function(z) function(b) 1)
    priors <- append(priors, function(b) {1/b})
  }
  names(priors) <- var.names
  
  data <- mutate(data, const.=rep(1,nrow(data))) #necessary for next step
  
  likelihood <- function(inputs, params) {
    k <- length(params) - 1
    params <- as.numeric(params)
    nv <- params[k+1]
    deviates <- apply(inputs, 1, function(pnt) 
      {pnt[1] - params[1:k] %*% pnt[2:(k+1)]})
    #print(summary(deviates))
    #print(nv)
    LLH <- sum(log(dt(deviates, nv)))
    LLH
  }
  
  #list of all parameters that show up in chain
  param.list <- as_tibble(rbind(starting)) 
  #let's chain
  accepts <- 0
  
  for (i in seq(iter)) {
    # calculate prior densities
    old.params <- param.list[i,]
    #print(old.params)
    proposed.params <- sapply(seq(k+1), function(x) {proposals[[x]](old.params[x])})
    names(proposed.params) <- var.names
    #print(proposed.params)
    old.prior.prob <- 1
    new.prior.prob <- 1
    for (j in seq(k)) {
      #old.prior.prob <- priors[[j]](old.params[j]) #prior prob of parameters, NOT USED
      #new.prior.prob <- priors[[j]](proposed.params[j])
    }
    alpha.numer <- likelihood(data, proposed.params)
    alpha.denom <- likelihood(data, old.params) #can be optimized to use LLH from last iteration
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
      new.row <- as_tibble(proposed.params)
      accepts <- accepts + 1
    } else {
      new.row <- old.params
    }
    #print(new.row)
    param.list <- param.list %>% add_row(new.row)
    #print("Full list:")
    #print(param.list)
    if (i %% 50 == 0) {
      cat(paste0("\r", "Number of Accepts / Proposals: ", accepts, " / ", i))
    }
  }
  cat("\nFinished")
  param.list[(burn+1):(iter+1),]
}


auto_mc <- function(data) {
  k <- ncol(data) #number of variables, not including intercept
  print("First running a test chain of 2000 samples (burning first 500):")
  test.chain <- metro_mc(data, iter=2000, burn=50)
  summary(test.chain)
  test.accepts <- nrow(unique(test.chain[,1]))
  test.accepts
  test.means <- apply(test.chain, 2, mean)
  test.var <- sapply(test.means, function(t) max(abs(t)/1.6, 1))
  if (test.accepts >= 3) {
    test.var <- apply(test.chain, 2, function(m) {(max(m)-min(m))/sqrt(3)})
    # i know this is messy, but using ifelse wasn't working, will try fixing it later
  }
  print(test.means)
  print(test.var)
  
  zone_in <- function(data.tib, last.chain, means, varis, batches=5, batch.size=1000, batch.burn=50) {
    if (batches == 0) {
      return(last.chain)
    } else {
      chain.proposals <- lapply(seq(k), function(z) function(old) old + rnorm(1, 0, varis[z]))
      chain.proposals <- append(chain.proposals, function(old) {1+(old-1)*exp(rcauchy(1, 0, varis[k+1]/5))})
      
      new.chain <- metro_mc(data, iter=batch.size, burn=batch.burn, proposals=chain.proposals, starting=means)
      
      chain.accepts <- nrow(unique(new.chain[,1]))
      chain.means <- apply(new.chain, 2, mean)
      chain.varis <- sapply(chain.means, function(t) max(abs(t), 1))
      if (chain.accepts >= 3) {
        chain.var <- apply(new.chain, 2, function(m) {(max(m)-min(m))/sqrt(3)})
        # i know this is messy, but using ifelse() wasn't working, will try fixing it later
      }
      print(rbind(chain.means, chain.varis))
      return(zone_in(data.tib, new.chain, chain.means, chain.varis, batches-1, batch.size, batch.burn))
    }
  }
  
  batching.test <- zone_in(data, test.chain, test.means, test.var, batch.size = 2000)
  
  chain1.proposals <- lapply(seq(k), function(z) function(old) old + rnorm(1, 0, test.var[z]))
  chain1.proposals <- append(chain1.proposals, function(old) {1+(old-1)*exp(rcauchy(1, 0, test.var[k+1]/5))})
  print("Running first chain with updated proposal distributions")
  chain1 <- metro_mc(data, iter=1000, burn=50, proposals=chain1.proposals, starting=test.means)
  summary(chain1)
  
  chain1.accepts <- nrow(unique(chain1[,1]))
  chain1.means <- apply(chain1, 2, mean)
  chain1.var <- sapply(chain1.means, function(t) max(abs(t), 1))
  if (test.accepts >= 3) {
    chain1.var <- apply(chain1, 2, function(m) {(max(m)-min(m))/sqrt(3)})
    # i know this is messy, but using ifelse wasn't working, will try fixing it later
  }
  chain1.means
  chain1.var
  
  chain2.proposals <- lapply(seq(k), function(z) function(old) old + rnorm(1, 0, chain1.var[z]/2))
  chain2.proposals <- append(chain2.proposals, function(old) {1 + (old-1)*exp(rcauchy(1, 0, chain1.var[k+1]))})
  print("Running first chain with updated proposal distributions")
  chain2 <- metro_mc(data, iter=1000, burn=50, proposals=chain2.proposals, starting=chain1.means)
  summary(chain2)
}

auto_mc(data)
