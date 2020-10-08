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
  
  if (is.null(starting)) {starting <- c(rep(0,k-1), mean(unlist(data[,1])), 1)}
  names(starting) <- c(names(data[,2:k]), "intercept", "nv")
  if (proposals == "rwalk") {
    proposals <- sapply(seq(k), function(t) function(old) {old+rnorm(1,0,10)})
    proposals <- append(proposals, function(old) {old*exp(rnorm(1,0,5))}) # last one is for the precision, tau
  }
  names(proposals) <- names(starting)
  if (priors == "vague") {
    priors <- sapply(seq(k), function(z) function(b) 1)
    priors <- append(priors, function(b) {1/b})
  }
  names(priors) <- names(starting)
  
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
    if (i %% 20 == 0) {
      cat(paste0("\r", "Number of Accepts / Proposals: ", accepts, " / ", i))
    }
  }
  param.list[(burn+1):iter,]
}


auto_mc <- function(data) {
  test.chain <- metro_mc(data, iter=2000, burn=500)
  test.accepts <- nrow(unique(test.chain[,1]))
  test.accepts
  test.means <- apply(test.chain, 2, mean)
  test.var <- sapply(test.means, function(t) max(t,1))
  if (test.accepts >= 20) {
    test.var <- apply(test.chain, 2, var)
    # i know this is messy, but using ifelse wasn't working, will try fixing it later
  }
  
  
}
