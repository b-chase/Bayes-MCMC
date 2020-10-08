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
    proposals <- append(proposals, function(t) rexp(1,0.05)) # last one is for the precision, tau
  }
  names(proposals) <- names(starting)
  if (priors == "vague") {
    priors <- sapply(seq(k+1), function(z) function(b) 1)
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
    proposed.params <- as.numeric(sapply(seq(k+1), function(x) {proposals[[x]](param.list[i,x])}))
    prior.prob <- 1
    for (j in seq(k)) {
      prior.prob <- prior.prob * priors[[j]](param.list[i,j]) #prior prob of params
    }
    alpha.numer <- likelihood(data, proposed.params)
    alpha.denom <- likelihood(data, param.list[i,])
    #print(alpha.numer)
    #print(alpha.denom)
    alpha <- alpha.numer - alpha.denom
    #print(alpha)
    if (log(runif(1)) < alpha) {
      param.list <- rbind(param.list, proposed.params)
      accepts <- accepts + 1
    } else {
      param.list <- rbind(param.list, param.list[i,])
    }
    if (i %% 20 == 0) {
      cat(paste0("\r", "Number of Accepts / Proposals: ", accepts, " / ", i))
    }
  }
  param.list[(burn+1):iter,]
}

test.chain <- metro_mc(data, iter=12000, burn=2000)
par(mfcol=c(2,2))
traceplot(mcmc(test.chain))


