data {
  int y_sample;
  int n_sample;
  real mu_spec;
  real mu_sens;
  real<lower = 0> sigma_spec;
  real<lower = 0> sigma_sens;
}
parameters {
  real<lower = 0, upper = 1> p;
  real<lower = 0, upper = 1> spec;
  real<lower = 0, upper = 1> sens;
}
model {
  real p_sample = p * sens + (1 - p) * (1 - spec);
  y_sample ~ binomial(n_sample, p_sample);
  spec ~ normal(mu_spec, sigma_spec);
  sens ~ normal(mu_sens, sigma_sens);
}
