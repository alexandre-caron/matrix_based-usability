
# Function estimates ------------------------------------------------------

fn_estimates <- function(d, screen_limit=50){
  n <- nrow(d)
  ms <- colSums(d)[colSums(d) > 0]
  j <- sum(ms>0)
  p <- sum(ms)/(n*j)
  sg <- sum(ms %in% 1)
  naive <- fit.naive(n=n, j=j, p=p)
  gt <- fit.gt(n=n, j=j, p=p, singletons=sg)
  dd <- fit.dd(n=n, j=j, p=p, singletons=sg)
  lnbzt <- tryCatch(expr = fit.LNBzt(ms = ms, n = n), error=function(e){list(zt=NA, mu=NA, sd=NA)})
  bayes <- tryCatch(expr = heterogeneous_bayes(j=j, n=n, ms=ms, M=ncol(d)+screen_limit), error = function(e){NA})
  return(tryCatch(
    expr = c(m_naive=naive, m_gt=gt[1], m_dd=dd[1], m_lnbzt=dtotal.LNBfit(lnbzt), m_bayes=bayes,
             p=p, p_gt=gt[2], p_dd=dd[2], mu_lnbzt=lnbzt$mu, sigma_lnbzt=lnbzt$sd), 
    error = function(e){c(m_naive=NA, m_gt=NA, m_dd=NA, m_lnbzt=NA, m_bayes=NA,
                          p=NA, p_gt=NA, p_dd=NA, mu_lnbzt=NA, sigma_lnbzt=NA)
      NA}
  ))
}


# Naive -------------------------------------------------------------------

fit.naive <- function(n, j, p){
  m <- j/(1-(1-p)^n)  
  return(unname(m))
}


# Good-Turing -------------------------------------------------------------

fit.gt <- function(n, j, p, singletons){
  p_gt <- p/(1+singletons/j)
  m_gt <- j/(1-(1-(p/(1+singletons/j)))^n)
  return(unname(c(m_gt,p_gt)))
}


# Double Deflation --------------------------------------------------------

fit.dd <- function(n, j, p, singletons){
  p_lewis <- 1/2*(p/(1+singletons/j))+1/2*((p-1/n)*(1-1/n))
  m_lewis <-  j/(1-(1-p_lewis)^n)
  return(unname(c(m_lewis,p_lewis)))
}


# Logit Normal Binomial zero truncated (LNBzt) Schmettow ------------------

## Logit, Inverse Logit, Density & Random
logit<-function(x) log(x/(1-x))
logitinv <- function(x) 1/(1+exp(-x))
dlogitnorm<-function(x,mu,sigma) dnorm(logit(x),mu,sigma)/(x*(1-x))
rlogitnorm<-function(r,mu,sigma) plogis(rnorm(r,mu,sigma))


rlnbinom<-function(r,size,m,s) rbinom(r,size,rlogitnorm(r,m,s))

## Helper Function
freq<-function(ms, occ = NULL){
  if(is.null(occ)) occ = unique(ms)
  sapply(occ, function(x) sum(ms == x))
}

## LNB probability distribution function
dlnbinom.mc <- function (x, size, m, s, nruns = 1e+05){
  f = function(x, size, m, s) freq(rlnbinom(nruns, size, m, s), x)/nruns
  return(mapply(f, x, size, m, s))
}

dlnbinom <- function(x, size, m, s){
  x = as.integer(x)
  size = as.integer(size)
  f = function(P, x, size, m, s) (1 - P)^(size - x - 1) * P^(x - 1) * exp(-(logit(P) - m)^2/(2 * s^2))
  F = function(x, size, m, s) {
    if (s > 0) return(integrate(f, 0, 1, x, size, m, s, rel.tol = .Machine$double.eps^0.25, subdivisions = 100)$value * choose(size, x)/(sqrt(2 * pi * s^2)))
    if (s == 0) return(dbinom(x, size, exp(m)))
    if (s < 0) return(0)
  }
  tryCatch(mapply(F, x, size, m, s), error = function(cond) {return(dlnbinom.mc(x, size, m, s, nruns = 10000))})
}

## LNBzt probability distribution function
dlnbinom.zt<-function(x,size,m,s){
  d0<-(1-dlnbinom(0,size,m,s))
  (x>0)*dlnbinom(x,size,m,s)/d0
}


## (Negative) LogLikelihood LNBzt
nloglik.lnbinom.zt<-function(p,n,K){ 
  K<-K[K>0]
  if(p[2]>=0){
    range<-c(c(1:n))
    frq<-freq(K,range)
    filter<-frq>0
    ml<-log(dlnbinom.zt(range[filter],n,p[1],p[2]))
    nloglik<-(-sum(ml*frq[filter]))
    return(nloglik)
  }else{return(-Inf)}
}


## Optim
fit.LNBzt<-function(ms,n,startval=c(-1,2)){
  ms = unlist(ms)
  maxLik<-optim(fn=nloglik.lnbinom.zt, par=startval,K=ms, n=n, method="BFGS")
  LNBfit<-list(nlogLik=maxLik$value, mu=maxLik$par[1], sd=maxLik$par[2],n=n, discovered=sum(ms>0), ms=ms)
  LNBfit$AIC=2*(-LNBfit$nlogLik)+2*2
  LNBfit$zt = T
  LNBfit <- structure(LNBfit, class = c("LNBfit"))
  return(LNBfit)
}

dhat.LNBfit <- function(object) ifelse(object$zt, 1-dlnbinom(0,object$n,object$mu, object$sd), sum(object$ms>0)/length(object$ms))
dtotal.LNBfit <- function(object) ifelse(object$zt, object$discovered / dhat(object), length(object$ms))
dnull.LNBfit <- function(object) ifelse(object$zt, object$discovered/dhat(object) - object$discovered, sum(object$ms==0))
dhat <- function(x) UseMethod("dhat",x)
dnull <- function(x) UseMethod("dnull",x)
dtotal <- function(x) UseMethod("dtotal",x)



# Bayesian ----------------------------------------------------------------

heterogeneous_bayes <- function(j,n,ms, M, h = 5, mu_prior_mean = 0, mu_prior_sd = 10000, 
                                 s2_prior_a = 0.5, s2_prior_b = 0.5) {
  simu <- NULL
  lp <- rep(-Inf, M)
  for (k in j:M){
    x <- c(ms, rep(0, k - j)) 
    usability_dat <- list(n = n, m = k, x = x, mu_prior_mean = mu_prior_mean, 
                          mu_prior_sd = mu_prior_sd, s2_prior_a = s2_prior_a,
                          s2_prior_b = s2_prior_b)
    fit <- sampling(object = model,  data = usability_dat, refresh = 0)
    lp[k] <- bridge_sampler(fit,silent = TRUE)$logml
    lp[k] <- lp[k] + lchoose(k,j)
    if ((k >= (h+j)) & (lp[k - h] > max(lp[(k - h + 1):k]))){
      k = k - h
      break
    } 
  }
  return(k)
}


heterogeneous_bayes_ci <- function(j,n,ms, M, h = 5, mu_prior_mean = 0, mu_prior_sd = 10000, 
                                   s2_prior_a = 0.5, s2_prior_b = 0.5) {
  
  simu <- NULL
  lp <- rep(-Inf, M)
  
  init <- 0 
  for (k in j:M){
    x <- c(ms, rep(0, k - j))
    usability_dat <- list(n = n, m = k, x = x, mu_prior_mean = mu_prior_mean, 
                          mu_prior_sd = mu_prior_sd, s2_prior_a = s2_prior_a,
                          s2_prior_b = s2_prior_b)
    fit <- sampling(object = model,  data = usability_dat, refresh = 0)
    lp[k] <- bridge_sampler(fit,silent = TRUE)$logml
    lp[k] <- lp[k] + lchoose(k,j)
    mu = do.call("c", lapply(fit@sim$samples, function(x) x$mu))
    s2 = do.call("c", lapply(fit@sim$samples, function(x) x$s2))
    simu = rbind(simu, cbind(mu = mu, s2 = s2, m = k))
  }
  
  lp <- lp - max(lp)
  posterior_m <- prop.table(exp(lp)) 
  posterior_m <- cbind(m = 1:M, posterior = posterior_m)
  return(list(m = which.max(posterior_m[,2]), posterior_m = posterior_m, simu = simu))
}


# Confidence Interval -----------------------------------------------------

fn_sim_homogeneous <- function(m, p, n){
  d <- matrix(rbinom(n=n*m, size=1, prob=p), nrow=n, ncol=m)
  d <- as.matrix(d[,colSums(d)>0])
  if(ncol(d)>0){return(d)}else{return(NULL)}
}

fn_sim_heterogeneous <- function(m, mu, sigma, n){
  pl <- logitinv(rnorm(n = m, mean = mu, sd = sigma))
  d <- sapply(pl, function(x) rbinom(n = n, size = 1, prob = x))
  d <- as.matrix(d[,colSums(d)>0])
  if(ncol(d)>0){return(d)}else{return(NULL)}
  return(d)
}

boot_naive <- function(m, p, n, n_boot = 1000){
  sapply(X = 1:n_boot, FUN = function(x){
    d <- fn_sim_homogeneous(m=m, p=p, n=n)
    if(is.null(d)){
      return(NA)
    }else{
      return(fit.naive(n=nrow(d), j=ncol(d), p = mean(d)))
    }
  })
}

boot_gt <- function(m, p, n, n_boot = 1000){
  sapply(X = 1:n_boot, FUN = function(x){
    d <- fn_sim_homogeneous(m=m, p=p, n=n)
    if(is.null(d)){
      return(NA)
    }else{
      return(fit.gt(n=nrow(d), j=ncol(d), p = mean(d), singletons = sum(colSums(d) %in% 1))[1])
    }
  })
}

boot_dd <- function(m, p, n, n_boot = 1000){
  sapply(X = 1:n_boot, FUN = function(x){
    d <- fn_sim_homogeneous(m=m, p=p, n=n)
    if(is.null(d)){
      return(NA)
    }else{
      return(fit.dd(n=nrow(d), j=ncol(d), p = mean(d), singletons = sum(colSums(d) %in% 1))[1])
    }
  })
}

boot_lnbzt <- function(m, mu, sigma, n, n_boot = 1000){
  sapply(X = 1:n_boot, FUN = function(x){
    if(is.na(m)|is.na(mu)|is.na(sigma)){
      d <- NULL
    }else{
      d <- fn_sim_heterogeneous(m=m, mu=mu, sigma=sigma, n=n)
      }
    if(is.null(d)){
      return(NA)
    }else{
      return(tryCatch(expr = dtotal.LNBfit(fit.LNBzt(ms=colSums(d), n=n)), error=function(e){NA}))
    }
  })
}



