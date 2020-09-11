data {
  int<lower=0> M;
  vector[2] y[M];
}

parameters {
  corr_matrix[2] Omega;
  vector[2] mu;
  vector<lower=0>[2] sigma;
}
transformed parameters {
  cov_matrix[2] Sigma;
  Sigma = quad_form_diag(Omega, sigma);
}
model {
  sigma ~ cauchy(0, 5);
  Omega ~ lkj_corr(1);
  mu ~ normal(0,5);
  y ~ multi_normal(mu, Sigma);}


