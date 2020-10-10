require(tidyverse)
library(coda)

n <- 500
x <- runif(n, -10, 30)
z <- runif(n, -30, -12)
b1 <- 42
b1.alt <- -24
b2 <- -0.02
a <- -13

prob.fail <- 0
y <- ifelse(runif(n)>prob.fail,b1*x,b1.alt*x) + b2*z + a + rnorm(n, 0, 25)

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
  
  var.names <- c(names(data[,2:k]), "intercept", "tau")
  # k+1 terms if we count tau (precision)
  
  if (is.null(starting)) {starting <- c(rep(0,k), 0.000001)}
  names(starting) <- var.names
  if (is.character(proposals)) {
    proposals <- function(old.props) {c(old.props[1:k] + runif(k,-10,10), rgamma(1, 1, 20)+10^-100)}
  }
  if (is.character(priors)) {
    priors <- function(params) {c(rep(1,k), 1/params[k+1])}
  }
  
  data <- mutate(data, const.=rep(1,nrow(data))) #necessary for next step
  
  likelihood <- function(inputs, params, k) {
    tau <- 1/sqrt(params[k+1])
    deviates <- apply(inputs, 1, function(pnt) 
      {pnt[1] - params[1:k] %*% pnt[2:(k+1)]})
    #print(summary(deviates))
    #print(tau)
    LLH <- sum(log(dnorm(deviates, 0, tau)))
    LLH
  }
  
  test <- c(1.7, -0.6815, -69.74, 0.001)
  an <- likelihood(data, test, k)
  an
  started <- c(13.453, -10.0885, -379.4, 2.267*10^-5)
  ad <- likelihood(data, started, k)
  ad
  proposals(test)
  priors(test)
  exp(an-ad)
  
  #list of all parameters that show up in chain
  param.list <- as.matrix(rbind(starting)) 
  old.LLH <- likelihood(data, starting, k)
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
    new.LLH <- likelihood(data, proposed.params, k)
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
    if (i %% 500 == 0) {
      cat(paste0("\rNumber of Accepts / Proposals: ", accepts, " / ", i))
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
  as_tibble(param.list[(burn+2):(iter+1),])
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
  test.chain <- metro_mc(test.data, iter=2000, burn=500, starting=lm.guesses)
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
      
      new.chain <- metro_mc(data, iter=batch.size, burn=batch.burn, proposals=chain.proposals, starting=means)
      
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

my.chain <- auto_mc(data, batch.sizes = 20000, total.batches = 4)
densplot(mcmc(my.chain))
densplot
par(mfcol=c(2,2))
apply(my.chain, 2, hist)

first <- as.matrix(rbind(c(0.0001,1,0)))
tests <- first

'
for (i in seq(5000)) {
  h1 <- tests[i,1]
  h2 <- tests[i,2]
  test.prop <- function(old.params) {c(old.params[1:k] + rcauchy(k,0), 
                                       rgamma(1,h1, h2)+10^-100)}
  chain.i <- metro_mc(data, iter=2000, burn=100, starting = c(30, -0.1, 0, 0.01), proposals = test.prop)
  tests[i,3] <- nrow(unique(chain.i[,1]))
  tests <- rbind(tests, c(10^runif(1,-2,0.5), 10^runif(1,-5,3), 0))
  print(paste0("Latest test result! (", i, ")"))
  print(apply(chain.i,2,mean))
  print(signif(tests[i,]),4)
}
tib.tests <- as_tibble(tests)
names(tib.tests) <- c("h1", "h2", "accepts")

tib.tests %>% arrange(desc(accepts)) %>% mutate("e(x)"=h1*h2, "v(x)"=h1*h2^2)
ggplot(filter(tib.tests, accepts>10)) +
  geom_point(mapping=aes(h1, h2, size=accepts), position = "jitter")

filter(tib.tests, accepts>18) %>% arrange(desc(accepts))
'