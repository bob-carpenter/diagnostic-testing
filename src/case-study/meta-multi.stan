data {
  int<lower = 0> K;                     // number of test sites
  int<lower = 1, upper = 3> status[K];  // status of subjects at site k
  int pos_tests[K];                     // number positive tests at site k
  int tests[K];                         // number of tests at site k
}
transformed data {
  int NEG = 1;
  int POS = 2;
  int UNK = 3;
}
parameters {
  matrix[2, 2] L_Omeg;
  // cholesky_corr[2] L_Omega;             // corr of sens and spec
  real<lower = 0> sigma[2];             // scale of variation in sens and spec
  vector[2] mu;                         // mean sens and spec
  vector<offset = mu,
         multiplier = diag_pre_multiply(sigma, L_Omega)>[2]
      logit_sens_spec[K];               // sens and spec at site k
  real<lower = 0, upper = 1> prev;      // prevalence
}
transformed parameters {
  vector[K] sens = to_vector(inv_logit(logit_sens_spec[ , 1]));
  vector[K] spec = to_vector(inv_logit(logit_sens_spec[ , 2]));
  vector[K] true_pos_rate = prev * sens;
  vector[K] false_pos_rate = (1 - prev) * (1 - spec);
  vector[K] pos_rate = true_pos_rate + false_pos_rate;
}
model {
  // hyperpriors
  mu ~ normal(2, 1);
  sigma ~ normal(0, { 0.2, 0.5 });
  L_Omega ~ lkj_corr_cholesky(4);

  // hierarchical prior
  logit_sens_spec ~ multi_normal_cholesky(mu, diag_pre_multiply(sigma, L_Omega));

  // prior
  prev ~ uniform(0, 1);

  // likelihood
  for (k in 1:K) {
    if (status[k] == NEG)
      pos_tests[k] ~ binomial(tests[k], 1 - spec[k]);
    else if (status[k] == POS)
      pos_tests[k] ~ binomial(tests[k], sens[k]);
    else if (status[k] == UNK)
      pos_tests[k] ~ binomial(tests[k], pos_rate[k]);
  }
}
generated quantities {
  vector[K] pos_pred_val = true_pos_rate / pos_rate;
}
