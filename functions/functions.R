# Required packages
library(rstan)
library(bridgesampling)

model <- stan_model('draw_mu_s2.stan') # compile the Stan program

# Input :
# - d : a discovery matrix
# - M : the maximal number of problem
# - h : the number of iterations to break after if no improvement of the 
# likelihood has been produced (default = 5)
# - ... : some optional inputs as the hyperparameters
# Ouput :
# - return the maximum a posteriori estimator of m (default), or posterior 
# distribution of m and resulting posterior sample of p(mu, s2|d, m) if 
# full_output = TRUE 
heterogeneous_bayes <- function(d, M, h = 5, mu_prior_mean = 0, mu_prior_sd = 10000, 
                                s2_prior_a = 0.5, s2_prior_b = 0.5, full_output = FALSE) {
  j = ncol(d)
  n = nrow(d)
  nl = colSums(d)
  if (full_output){
    simu = NULL
  }
  lp = rep(-Inf, M)
  for (k in j:M){
    # print(paste("k = ",k))
    x <- c(nl, rep(0, k - j)) # pad the vector with 0 for undiscovered problems
    usability_dat <- list(n = n, m = k, x = x, mu_prior_mean = mu_prior_mean, 
                          mu_prior_sd = mu_prior_sd, s2_prior_a = s2_prior_a,
                          s2_prior_b = s2_prior_b)
    fit <- sampling(object = model,  data = usability_dat, refresh = 0)
    # fit <- stan(file = 'draw_mu_s2.stan', data = usability_dat, 
    #            verbose = FALSE, open_progress = 0)
    lp[k] <- bridge_sampler(fit,silent = TRUE)$logml
    lp[k] <- lp[k] + lchoose(k,j)
    # print(paste0(k, "-",lp[k]))
    # It is eventually possible to get values of mu and s2 given m
    # by removing the 1000 first warmup iterations
    # But the resulting object can eventually be very big
    if (full_output){
      mu = do.call("c", lapply(fit@sim$samples, function(x) x$mu[1001:2000]))
      s2 = do.call("c", lapply(fit@sim$samples, function(x) x$s2[1001:2000]))   
      simu = rbind(simu, cbind(mu = mu, s2 = s2, m = k))
    } else {
      # If no amelioration during the h last iteration => stop
      if ((k >= (h+j)) & (lp[k - h] > max(lp[(k - h + 1):k]))){
        k = k - h
        break
      }      
    }
  }
  if (full_output){
    lp <- lp - max(lp)
    posterior_m <- prop.table(exp(lp))
    posterior_m <- cbind(m = 1:M, posterior = posterior_m)
    return(list(posterior_m = posterior_m, simu_mu_s2 = simu))
  } else {
    return(k)
  }
}
