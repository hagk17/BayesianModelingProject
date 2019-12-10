# MA 578, Bayesian Statistics
# Bayesian Linear Models
# y | beta, sigma2 ~ N(X * beta, sigma2 * I_n)
# beta ~ N(beta0, S_inv0 ^ (-1))
# sigma2 ~ Inv-scaled-Chisq(nu, tau2)

# [ Auxiliary ]
# sample from inverse (scaled) chi-square with parameters `nu` and `tau2`;
# nu_tau2 = nu * tau2 for convenience

#sample from an inverse chi square 
rinvchisq <- function (ns, nu, nu_tau2){ return(1 / rgamma(ns, nu / 2, nu_tau2 / 2)) }

#sets up the array to store the outcomes from the mcmc model 
#inputs: number of iterations, number of chains (default 1), parameter vector 
mcmc_array <- function (ns, nchains = 1, params) {
  nparams <- length(params) #number of parameters 
  #creates the 3d array to store the param values at each step of the chain
  array(dim = c(ns, nchains, nparams),
        dimnames = list(iterations = NULL,
                        chains = paste0("chain", 1:nchains),
                        parameters = params))
}


# [ Simple interface ]
bslm_fit3 <- function (y, x, prior_coef = NULL, prior_disp = NULL,
                      maxit = 25, epsilon = 1e-8) {
  nvars <- ncol(x) #number of variables 
  nobs <- nrow(x) #number of observations 
  dn <- colnames(x); #getting variable names 
  #if there are no names label x1,...,x_nvariables 
  if (is.null(dn)){ dn <- paste0("x", 1L:nvars) }
 # if (is.null(dn)){ dn <- paste0("x") }
  
  
  #if beta prior coefficients are not specified, 
  #set the mean and precision to be a list of zeros 
  if (is.null(prior_coef)){
    prior_coef <- list(mean = rep(0, nvars), precision = rep(0, nvars))
  }
  #if sigma^2 prior coefficients are not specified,
  #set degrees of freedom and scale parameters to 0 
  if (is.null(prior_disp)){
    prior_disp <- list(df = 0, scale = 0)
  }
  
  S_inv0 <- prior_coef$precision #set value of precision in beta prior 
  beta0 <- prior_coef$mean #set value of beta0 in beta prior
  
  #recalculate beta0 by multiplying by the precision for convenience 
  beta0 <- if (is.vector(S_inv0)) S_inv0 * beta0 else S_inv0 %*% beta0
  
  nu <- prior_disp$df #setting value of nu in sigma prior
  nu_tau2 <- nu * prior_disp$scale #multiplying nu and the scale together for convenience
  
  #sum((observed values-mean)^2) for intercept only 
  rss <- sum((y - mean(y)) ^ 2) 
  #calculates value of sigma^2 for expecation step 
  sigma2 <- (nu_tau2 + rss) / (nu + nobs)
  
  #This is the EM algorithm 
  #Goal: (beta, sigma2 = argmax{P(beta,sigma2 | y)})
  #Expecation Step: sigma2 from above 
  #Maximization Step: beta = (1/sigma2 * XTX+sigma0^{-1})^{-1}(1/sigma2 * XTy + beta0)
  # for the maximum number of iterations 
  for (iter in 1:maxit) {
    
    #setting up next beta = V_inv*z 
    #V_inv = (1/sigma2 * XTX+sigma0^{-1})^{-1}
    V_inv <- crossprod(x) / sigma2
    if (is.vector(S_inv0)){ # diagonal precision?
      diag(V_inv) <- diag(V_inv) + S_inv0
    }
    else{V_inv <- V_inv + S_inv0}
    #z = 1/sigma2 * XTy + beta0 
    z <- crossprod(x, y) / sigma2 + beta0
    
    #Steps for computing beta efficiently 
    #1) computing the cholesky factor of V_inv 
    C <- chol(V_inv)
    #2) compute z (which we already did)
    #3) compute Beta by using matrix magic 
    u <- backsolve(C, z, transpose = TRUE)
    coef <- drop(backsolve(C, u))
    
    #update the new rss and sigma2 (E step?)
    rss_new <- sum((y - drop(x %*% coef)) ^ 2)
    sigma2 <- (nu_tau2 + rss_new) / (nu + nobs)
    
    rel_error <- abs((rss_new - rss) / rss) 
    #checks that its not infinite and that it converged
    if (!is.infinite(rss_new) && (rel_error < epsilon)) # converged?
      break #break if did 
    rss <- rss_new # set the new rss and repeat 
  }
  names(coef) <- dn #name the coefficients the variable 
  list(coef = coef, sigma2 = sigma2, C = C) 
}


bslm_sample <- function (y, x, prior_coef = NULL, prior_disp = NULL,
                         chains = 4, iter = 2000, warmup = floor(iter / 2)) {
  nvars <- ncol(x) 
  nobs <- nrow(x)
  dn <- colnames(x)
  
  if (is.null(dn)) dn <- paste0("x", 1L:nvars)
  
  if (is.null(prior_coef))
    prior_coef <- list(mean = rep(0, nvars), precision = rep(0, nvars))
  
  if (is.null(prior_disp))
    prior_disp <- list(df = 0, scale = 0)
  
  S_inv0 <- prior_coef$precision
  
  beta0 <- prior_coef$mean
  
  beta0 <- if (is.vector(S_inv0)) S_inv0 * beta0 else S_inv0 %*% beta0
  
  nu <- prior_disp$df
  
  nu_tau2 <- nu * prior_disp$scale
  
  rss <- sum((y - mean(y)) ^ 2)
  sigma2 <- (nu_tau2 + rss) / (nu + nobs)
  sims <- mcmc_array(iter - warmup, chains, c(dn, "sigma2", "lp__"))
  for (chain in 1:chains) {
    for (it in 1:iter) {
      z <- crossprod(x, y) / sigma2 + beta0
      V_inv <- crossprod(x) / sigma2
      if (is.vector(S_inv0)) # diagonal precision?
        diag(V_inv) <- diag(V_inv) + S_inv0
      else
        V_inv <- V_inv + S_inv0
      C <- chol(V_inv)
      u <- backsolve(C, z, transpose = TRUE)
      coef <- drop(backsolve(C, u + rnorm(nvars)))
      rss <- sum((y - drop(x %*% coef)) ^ 2)
      sigma2 <- rinvchisq(1, nu + nobs, nu_tau2 + rss)
      lp <- -((nu + nobs) / 2 + 1) * log(sigma2) - .5 * (nu_tau2 + rss) / sigma2
      if (it > warmup)
        sims[it - warmup, chain, ] <- c(coef, sigma2, lp)
    }
  }
  sims
}

bslm_sample1 <- function (...) {
  s <- drop(bslm_sample(..., chains = 1, iter = 1)) # single sample
  head(s, -2) # remove sigma2 and lp__, leave coef only
}
