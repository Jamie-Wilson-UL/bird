// Weibull Survival Model for Bayesian Imputation 
// Based on Moghaddam et al. (2022) methodology

data {
  // Sizes
  int<lower=0> n_obs;                // number of observed events
  int<lower=0> n_cens;               // number right-censored

  // Data (time units as provided from R)
  vector<lower=0>[n_obs]  t_obs;
  vector<lower=0>[n_cens] t_cens;

  // Prior family: 1 = lognormal shape + Gamma scale (default), 2 = lognormal.
  int<lower=1, upper=2> prior_family;

  // Gamma prior for the positive Weibull scale parameter.
  real<lower=0> scale_prior_shape;
  real<lower=0> scale_prior_rate;

  // Lognormal priors: if log(shape) ~ Normal(mu_log_shape, sd_log_shape),
  // then shape ~ Lognormal(mu_log_shape, sd_log_shape), similarly for scale.
  real          mu_log_shape;
  real<lower=0> sd_log_shape;
  real          mu_log_scale;
  real<lower=0> sd_log_scale;

  // Weak upper bounds to prevent overflow in transformed params
  real<lower=1e-8> shape_upper;      
  real<lower=1e-8> scale_upper;      
}

parameters {
  // Directly sample on the positive scale with upper bounds 
  real<lower=1e-12, upper=shape_upper> shape;
  real<lower=1e-12, upper=scale_upper> scale;
}

model {
  // Priors
  if (prior_family == 1) {
    shape ~ lognormal(mu_log_shape, sd_log_shape);
    scale ~ gamma(scale_prior_shape, scale_prior_rate);
  } else {
    shape ~ lognormal(mu_log_shape, sd_log_shape);
    scale ~ lognormal(mu_log_scale, sd_log_scale);
  }

  // Likelihood
  target += weibull_lpdf(t_obs  | shape, scale);   // events: log f(t)
  target += weibull_lccdf(t_cens | shape, scale);  // censored: log S(t)
}

generated quantities {
  // Posterior imputations (T | T > t_cens) via inverse-CDF
  vector<lower=0>[n_cens] t_imputed;

  // Pointwise log-lik (for WAIC/LOO)
  vector[n_obs]  ll_obs;
  vector[n_cens] ll_cens;

  // Aggregates 
  real log_lik_obs   = weibull_lpdf(t_obs  | shape, scale);
  real log_lik_cens  = weibull_lccdf(t_cens | shape, scale);
  real log_lik_total = log_lik_obs + log_lik_cens;

  // Fill pointwise log-lik
  for (i in 1:n_obs)  ll_obs[i]  = weibull_lpdf(t_obs[i]  | shape, scale);
  for (i in 1:n_cens) ll_cens[i] = weibull_lccdf(t_cens[i] | shape, scale);

  // Impute: draw p ~ Uniform(F(tc), 1), then invert CDF
  for (i in 1:n_cens) {
    real u   = uniform_rng(0, 1);
    real Fc  = exp(weibull_lcdf(t_cens[i] | shape, scale));   // P(T <= t_cens)
    if (Fc < 1e-12)       Fc = 1e-12;
    if (Fc > 1 - 1e-12)   Fc = 1 - 1e-12;

    real p   = Fc + u * (1 - Fc);
    if (p < 1e-12)        p = 1e-12;
    if (p > 1 - 1e-12)    p = 1 - 1e-12;

    // Inverse Weibull: t = scale * (-log(1 - p))^(1/shape)
    real log1mp     = log1m(p);                 // stable for 1-p near 0
    real log_neglog = log(-log1mp);
    t_imputed[i]    = scale * exp(log_neglog / shape);

    // enforce right-censoring monotonicity
    if (t_imputed[i] <= t_cens[i])
      t_imputed[i] = t_cens[i] * (1 + 1e-12);
  }
}
