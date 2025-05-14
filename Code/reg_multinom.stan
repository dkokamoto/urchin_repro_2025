data {
  int<lower=1> N;  // number of observations
  array[N] int Y;  // response vector
  int N_pred;
  int K; // number of fixed effects parameters
  int G; // number of random effects parameters
  int<lower=2> ncat;  // number of categories
  matrix[N, K] X;  // fixed effects design matrix
  matrix[N, G] Z;  // random effects design matrix
  matrix[N_pred, K] X_pred;
  matrix[N_pred, G] Z_pred;
}

parameters {
  matrix<lower= 0> [K, ncat] b;  // fixed effects regression coefficients
  matrix [G, ncat - 1] z;  // random effects regression coefficients
  real<lower= 0> sigma;
}

transformed parameters {
  matrix[N, ncat] mu_logit;
  matrix[N, ncat] mu;
  matrix[K, ncat] beta;
  
  beta = log(b);
  mu_logit = X * beta +  Z * append_col(rep_vector(0.0, G),z);
 
  for(i in 1:N){
    mu[i] = to_row_vector(softmax(to_vector(mu_logit[i])));
  }
}

model {
  for(i in 1: K){
    target += gamma_lpdf(b[i] | 1, 1);
  }
  for (i in 1:G){
    target += std_normal_lpdf(z[i]);
  }
  target += student_t_lpdf(sigma | 3, 0, 0.1) - student_t_lccdf(0 | 3,0,0.1);
  for(i in 1:N){
    target += categorical_lpmf(Y[i] | to_vector(mu[i]));
  }
}

generated quantities {
  matrix[N_pred,ncat] mu_pred;
  vector[N] log_lik;
  
  for(i in 1:N){
    log_lik[i] = categorical_lpmf(Y[i] | to_vector(mu[i]));
  }
  
  {
    matrix[N_pred,ncat] mu_pred_logit = X_pred * beta +  Z_pred * append_col(rep_vector(0.0, G),z);
    for(i in 1:N_pred){
      mu_pred[i] = to_row_vector(softmax(to_vector(mu_pred_logit[i])));
    }
  }
}

