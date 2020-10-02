require(tidyverse)

linear_model <- function(yvals, xvals, b, b0, t) {
  #print("likelihood")
  #print(as.matrix(b))
  e <- unlist(yvals - as.matrix(xvals) %*% as.matrix(b) - b0)
  #prob <- prod(unlist(sapply(e,dnorm, 0, 1/sqrt(t))))
  prob <- prod(dnorm(unlist(e),0,1/sqrt(t)))
  if (is.nan(prob) || is.na(prob)) {
    print("Liklihood is NaN")
    print(cbind(yvals, xvals, e,dnorm(e,0,1/sqrt(t))))
    print(cbind(b, b0, t))
  }
  return(prob)
}

log_linear_model <- function(yvals, xvals, b, b0, t) {
  e <- unlist(log(yvals) - as.matrix(xvals) %*% as.matrix(b) - b0)
  #prob <- prod(unlist(sapply(e,dnorm, 0, 1/sqrt(t))))
  prob <- prod(dnorm(unlist(e),0,1/sqrt(t)))
  if (is.nan(prob) || is.na(prob)) {
    print("Liklihood is NaN")
    print(cbind(yvals, xvals, e,dnorm(e,0,1/sqrt(t))))
    print(cbind(b, b0, t))
  }
  return(prob)
}

"Metro_Hast" <- function(iter, dataset, betaprior=function(b){dnorm(b,0,100)}, 
                         beta0prior=function(b0){dnorm(b0,0,100)}, 
                         tauprior=function(t){1/t}, likelihood = linear_model, 
                         betaguess = 0, beta0guess = 0, tauguess = 0.001, 
                         betaprop = function(oldbeta,k){oldbeta + rnorm(k)}, 
                         beta0prop = function(oldbeta0){oldbeta0+rnorm(1)}, 
                         tauprop = function(oldtau){rgamma(1,1)}, 
                         return_alphas = FALSE) {
  #new version of metropolis hastings algorithm
  ydat <- dataset[,1]
  k = dim(as.matrix(dataset))[2]-1 #number of parameters
  xdat <- dataset[,2:(k+1)]
  
  #building outputs of the algorithm
  betas <- matrix(nrow=iter-1, ncol=k)
  betas <- rbind(betaguess,betas)
  #betas <- matrix(rep(ifelse(length(betaguess)==k, betaguess, rep(0,k)), iter), ncol = k)
  beta0s <- rep(beta0guess,iter)
  taus <- rep(tauguess,iter)
  
  #debug
  #print(betas)
  
  alphas <- rep(NA,iter)
  
  for (i in 2:iter) {
    if (i%%200==0) {cat(paste0("\r", i))}
    
    betas[i,] <- betaprop(betas[i-1,], k)
    beta0s[i] <- beta0prop(beta0s[i-1])
    taus[i] <- tauprop(taus[i-1])
    
    #rejection step
    u <- runif(1)
    
    #just metro right now, need to refix with metro hastings
    alpha_numer <- likelihood(ydat, xdat, betas[i,], beta0s[i], taus[i])*prod(betaprior(betas[i,]))*beta0prior(beta0s[i])*tauprior(taus[i])
    alpha_denom <- likelihood(ydat, xdat, betas[i-1,], beta0s[i-1], taus[i-1])*prod(betaprior(betas[i-1,]))*beta0prior(beta0s[i-1])*tauprior(taus[i-1])
    
    alpha <- alpha_numer / alpha_denom
    #print(paste0(i,": Numer: ", alpha_numer, " Denom: ", alpha_denom))
    alphas[i] <- alpha
    
    if (is.nan(alpha_numer)) {
      print("Numerator is NaN")
      return(tibble(betas, beta0s, taus, alphas))
    }
    if (is.na(alpha_numer)) {
      print("Numerator is NA")
      return(tibble(betas, beta0s, taus, alphas))
    }
    if (alpha_numer > alpha_denom & alpha_denom == 0) {
      print(paste0("Alpha denominator is 0 in iteration: ", i,"."))
      print(paste0("Numerator is ", alpha_numer))
      alpha <- 1
    }
    #if there is an error, break loop and return message
    if(is.na(alpha)) {
      print(paste0("Alpha is NA! On iteration ", str(i)))
      print(paste0("Numerator: ", alpha_numer, " Denominator: ", alpha_denom))
      return(tibble(betas, beta0s, taus, alphas))}
    #rejection step
    if (u > alpha) {
      betas[i,] <- betas[i-1,]
      beta0s[i] <- beta0s[i-1]
      taus[i] <- taus[i-1]
    }
  }
  
  postsample <- as_tibble(betas)
  postsample$beta0 <- beta0s
  postsample$tau <- taus
  if (return_alphas) {postsample$alpha <- alphas}
  
  return(postsample)
}

n <- 2000
x <- runif(n, -10, 30)
b <- 27
a <- -2
y <- b*x + a + rnorm(n, 0, 2)

lm(y~x, tibble(y,x))
