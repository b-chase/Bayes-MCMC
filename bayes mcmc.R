require(tidyverse)
library(coda)

lin_LLH <- function(inputs, params, k) {
  tau <- 1/sqrt(params[k+1])
  deviates <- apply(inputs, 1, function(pnt) 
  {pnt[1] - params[1:k] %*% pnt[2:(k+1)]})
  #print(summary(deviates))
  #print(tau)
  LLH <- sum(log(dnorm(deviates, 0, tau)))
  LLH
}

metro_mc <- function(data, iter = 10000, burn = 500, starting=NULL, 
                     proposals="rwalk", priors="vague", 
                     add.ones=TRUE, likelihood=lin.LLH) {
  # data is a table with dependent variable in first column
  # mmc_spec is a list of two functions, priors and LH
  # prior is a function for drawing potential values of the parameters
  # LH is the likehood, a function that gives p density for a given error
  # iter is number of MC steps to take
  # burn is burn-in, number of initial steps to ignore
  
  # next step - make sure priors and proposals return correct size vectors, or else return error
  # also make sure likelihood is a log likelihood, and returns correctly
  # also make sure starting params are workable and right size
  
  k <- ncol(data) #n-1 independent variables plus a constant
  
  var.names <- c(names(data[,2:k]), "exp", "intercept", "tau")
  # k+1 terms if we count tau (precision)
  
  if (is.null(starting)) {
    starting <- c(rep(0,k), 0.000001)
    var.names <- c(names(data[,2:k]), "intercept", "tau")
    names(starting) <- var.names
  } else {
    var.names <- names(starting)
    k <- length(starting)-1 # sometimes we won't have a normal linear model
  }
  if (is.character(proposals)) {
    proposals <- function(old.props) {c(old.props[1:k] + runif(k,-10,10), rgamma(1, 1, 20)+10^-100)}
  }
  if (is.character(priors)) {
    priors <- function(params) {c(rep(1,k), 1/params[k+1])}
  }
  if (add.ones) {data <- mutate(data, const.=rep(1,nrow(data)))}
  
  
  # debug
  #test <- c(1.7, -0.6815, -69.74, 0.001)
  #an <- likelihood(data, test, k)
  #an
  #started <- c(13.453, -10.0885, -379.4, 2.267*10^-5)
  #ad <- likelihood(data, started, k)
  #ad
  #proposals(test)
  #priors(test)
  #exp(an-ad)
  
  #list of all parameters that show up in chain
  param.list <- as.matrix(rbind(starting)) 
  old.LLH <- likelihood(data, starting, k)
  #let's chain
  accepts <- 0
  last.accepts <- 0
  likelihoods <- c(old.LLH)
  old.params <- starting

  for (i in seq(iter)) {
    # calculate prior densities
    #print(old.params)
    proposed.params <- proposals(old.params)
    #print(proposed.params)
    old.prior.prob <- prod(priors(old.params))
    new.prior.prob <- prod(priors(proposed.params))
    new.LLH <- likelihood(data, proposed.params, k)
    alpha.numer <- new.LLH + log(new.prior.prob)
    alpha.denom <- old.LLH + log(old.prior.prob)
    #print(alpha.numer)
    #print(alpha.denom)
    alpha <- alpha.numer - alpha.denom
    if (is.na(alpha)) {
      # output error information
      print(proposed.params)
      print(c("Prior, LLH, Numer", new.prior.prob, new.LLH, alpha.numer))
      print(old.params)
      print(c("Prior, LLH, Denom", old.prior.prob, old.LLH, alpha.denom))
    }
    #print(alpha)
    new.row <- c()
    if (log(runif(1)) < alpha) {
      new.row <- proposed.params
      accepts <- accepts + 1
      old.LLH <- new.LLH
      old.params <- proposed.params
    } else {
      new.row <- old.params
    }
    #print(new.row)
    param.list <- rbind(param.list, new.row)
    likelihoods <- c(likelihoods, old.LLH)
    #print("Full list:")
    #print(param.list)
    check.div <- 100
    if (i %% check.div == 0) {
      new.accepts <- accepts-last.accepts
      last.accepts <- accepts
      cat(paste0("\rNumber of Accepts / Proposals: ", accepts, " / ", i, 
                 "    -    Current Rate of Accepts: ", round(100*new.accepts/check.div, 1), "%"))
      
      #print("Last Accepted:")
      #print(new.row)
      #print(alpha.denom)
      #print("Last Proposed and LLH:")
      #print(proposed.params)
      #print(alpha.numer)
      #print("Last Alpha (difference)")
      #print(exp(alpha))
    }
  }
  cat("\nFinished\n")
  cat(paste0("There were ",accepts, " accepted proposals.\n"))
  names(param.list) <- var.names
  chain.out <- as_tibble(param.list[(burn+2):(iter+1),])
  names(chain.out) <- var.names
  list(chain=chain.out, llh.list=likelihoods[(burn+2):(iter+1)])
}


auto_mc <- function(data, batch.sizes=5000, total.batches=5) {
  k <- ncol(data) #number of variables, not including intercept
  
  var.names <- names(data)
  lin.model <- lm(eval(paste0(var.names[1],"~", paste0(var.names[2:k],collapse="+"))), data[1:20,])
  print("Linear model of first 20 terms:")
  print(summary(lin.model))
  lm.guesses <- c(lin.model$coefficients[c(2:k,1)], 1/var(lin.model$residuals))
  names(lm.guesses) = c(names(data[2:k]), "intercept", "tau")
  test.data <- data[21:nrow(data),]
  test.proposals <- function(old.params) {c(rnorm(k, lm.guesses[1:k]), 
                                            rgamma(1, 0.05, 0.02))}
  
  print("First running a test chain of 2000 samples (burning first 500):")
  test.chain <- metro_mc(test.data, iter=2000, burn=500, starting=lm.guesses)$chain
  summary(test.chain)
  test.accepts <- nrow(unique(test.chain[,1]))
  test.accepts
  test.means <- apply(test.chain, 2, mean)
  test.var <- sapply(test.means, function(t) max(abs(t)/5, 1))
  if (test.accepts >= 3) {
    test.var <- apply(test.chain, 2, function(m) {(max(m)-min(m))/sqrt(3)})
    # i know this is messy, but using ifelse wasn't working, will try fixing it later
  }
  print(as_tibble(signif(rbind(test.means, test.var)),5))
  Sys.sleep(5)
  
  zone_in <- function(data.tib, last.chain, means, varis, batches=5, batch.size=5000, batch.burn=500) {
    k <- ncol(data.tib)
    if (batches == 0) {
      return(last.chain)
    } else {
      cat(paste0("\nBatches remaining: ", batches, "\n"))
      chain.proposals <- function(old.params) {c(old.params[1:k] + rnorm(k,0,sqrt(varis[1:k])), 
                                                 old.params[k+1]*exp(rnorm(1,0,sqrt(varis[k+1]))))}
      
      new.chain <- metro_mc(data, iter=batch.size, burn=batch.burn, proposals=chain.proposals, starting=means)$chain
      
      chain.accepts <- nrow(unique(new.chain[,1]))
      chain.means <- apply(new.chain, 2, mean)
      chain.varis <- varis/1.6
      if (chain.accepts >= 5) {
        chain.varis <- apply(new.chain, 2, function(m) {(max(m)-min(m))/sqrt(3)})
      }
      print(as_tibble(signif(rbind(chain.means, chain.varis)),5))
      Sys.sleep(5)
      
      return(do.call(zone_in, args=list(data.tib, new.chain, chain.means, chain.varis, batches-1, batch.size, batch.burn)))
    }
  }
  
  last.chain <- zone_in(test.data, test.chain, test.means, test.var, batches=total.batches, batch.size=batch.sizes, batch.burn=0)
  par(mfcol=c(2,ceiling(k/2)+1))
  traceplot(mcmc(last.chain))  
  
  last.chain
}
