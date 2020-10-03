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

n <- 5000
x <- runif(n, -10, 30)
z <- runif(n, -30, -12)
b1 <- 27
b2 <- 0.04
a <- -2
y <- b1*x + b2*z + a + rnorm(n, 0, 1)

lm(y~x+z, tibble(y,x,z))

metromc_spec <- function(k, type=list('norm', 't', 'exp'), params=NA, 
                         intercept_prior=function(x){rnorm(x,0,)}
                         err_type=c('norm', 't', 'exp'), tau0=0.1) {
  # Specify priors and likelihood for doing metropolis algorithm MCMC
  # k is number of independent variables, not including intercept
  # params is a list of vectors containing parameters for each prior dist
  # (mean, sd) for norm, (df, ncp) for t, ()
  # tau0 is the initial value of the precision
  
  #defining priors
  priors <- list()
  for (i in seq(k+1)) {
    if (type[i]=='norm') {
      if (is.NA(params)) {
        mean <- 0
        sd <- 1
      }
      else {
        mean <- params[[i]][1]
        sd <- params[[i]][2]
      }
      priors <- append(priors, function(x) {rnorm(x, mean, sd)})
    }
    if (type[i]=='t') {
      if (is.NA(params)) {
        df <- 5
        ncp <- 0
      }
      else {
        df <- params[[i]][1]
        ncp <- params[[i]][2]
      }
      priors <- append(priors, function(x) {rt(x, df, ncp)})
    }
    if (type[i]=='exp') {
      if (is.NA(params)) {
        rate <- 1
      }
      else {
        rate <- params[[i]][1]
      }
      priors <- append(priors, function(x) {rexp(x, rate)})
    }
  }
  
  priors <- append(priors, intercept_prior) #intercept prior
  
  if (err_type=='norm') {
    priors <- append(priors, function(x) {rgamma(x, 1)}) #tau prior, eventually replace
    LHfunc <- function(x, t) {dnorm(x, 0, 1/t)}
  }
  return(list(priors, LHfunc))
}

auto_metromc_spec <- function(data) {
  # data is a table with dependent variables in first column
  # data doesn't contain column of '1's for intercept
  
  data <- tibble(y,x,z)
  k <- ncol(data) - 1 # number of independent variables
  n <- nrow(data) # number of datapoints
  red.n <- floor(n/4) #reduced n
  fixed.data <- data[1:red.n,] %>% mutate("intercept" = rep(1,red.n)) # adding row of ones for doing analysis
  
  priors <- list()
  new.prior <- function(x) {force(x); function(z) rnorm(z, x[1], sqrt(x[2]))}
  prior.means <- c()
  
  #getting summary stats for usage in error priors
  lean.y <- unlist(fixed.data[,1]) #for finding intercept later
  #y.mean <- mean(y)
  y.var <- var(lean.y)
  
  for (i in seq(k)) {
    paired.data <- fixed.data[,c(1,i+1)]
    lin.form <- paste0(names(paired.data)[1],"~",names(paired.data)[2])
    lin.m <- lm(lin.form, paired.data)
    beta.guess <- lin.m$coefficients[2]
    beta.g.var <- (sum(lin.m$residuals^2)/(red.n - (k+1))) / sum(paired.data[,2]^2)
    
    #print(beta.guess)
    #print(beta.g.var)
    prior.means <- c(prior.means, beta.guess)
    lean.y <- lean.y - (beta.guess * unlist(fixed.data[,i+1]))
    
    #need to fix this problem right here!!!
    
    
    priors <- append(priors, new.prior(c(beta.guess, 100*beta.g.var))) #go for vague priors
  }
  intercept.guess <- mean(lean.y)
  variance.guess <-  var(lean.y)
  df.guess <- 2*variance.guess/(variance.guess-1)
  
  intercept.prior <- new.prior(c(intercept.guess, 100*variance.guess))
  priors <- append(priors, intercept.prior)
  
  precision.prior <- function(n) rgamma(n, 1/(5*y.var), 5)
  priors <- append(priors, precision.prior)
  
  names(priors) <- c(paste0(""))
  priors
}

metromc <- function(data, mmc_spec, iter = 2000, burn = 500) {
  # data is a table with dependent variable in first column
  # mmc_spec is a list of two functions, priors and LH
    # prior is a function for drawing potential values of the parameters
    # LH is the likehood, a function that gives p density for a given error
  # iter is number of MC steps to take
  # burn is burn-in, number of initial steps to ignore
   
}

post_approx <- function() {}
