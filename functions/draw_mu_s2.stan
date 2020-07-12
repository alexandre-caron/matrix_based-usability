// Stan code for sampling mu and s2 given m and the data x[m]
// Sampling of m is still a problem ... 
data {
  int<lower=0> n;         // number of subjects 
  int<lower=0> m;         // number of problems 
  int<lower=0> x[m];      // number of dectection of each problem
  
  // hyperparameters
  real mu_prior_mean;
  real<lower=0> mu_prior_sd;
  real<lower=0> s2_prior_a;
  real<lower=0> s2_prior_b;
}
parameters {
  real mu;                // mean of the logit_p
  real<lower=0> s2;       // variance of logit_p
  real logit_p[m];        // vector of logit_p
}
model {
  target += normal_lpdf(mu | mu_prior_mean, mu_prior_sd);    // prior for mu
  target += inv_gamma_lpdf(s2 | s2_prior_a, s2_prior_b);  // prior for s2
  for (j in 1:m) {
    target += normal_lpdf(logit_p[j] | mu, sqrt(s2)); // logit_p given mu and s
    target += binomial_lpmf(x[j] | n, inv_logit(logit_p[j]));   // log-likelihood of the data
  }
}
