data {
  int<lower=0> M;
  vector[M] y;
  cov_matrix[2] SigmaEssential;
  vector[2] meanEssential;
  vector[2] meanNE;
  cov_matrix[2] SigmaNE;
}

parameters {
  vector[2] muSigma;
}
model {
  target += log_mix(0.1,multi_normal_lpdf(muSigma|meanEssential,SigmaEssential),
  multi_normal_lpdf(muSigma|meanNE,SigmaNE));
  y ~ normal(muSigma[1],muSigma[2]);}

