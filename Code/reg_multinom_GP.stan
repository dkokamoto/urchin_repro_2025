data {
  int<lower=1> N;  // number of observations
  array[N] int Y;  // response vector
  int N_pred;
  int K; // number of fixed effects parameters
  int G; // number of random effects parameters
  int<lower=2> ncat;  // number of categories
  matrix[N, K] X;  // fixed effects design matrix
  matrix[N, G] Z;  // random effects design matrix
  int Ngp;
  array[2] int N_gp_ID;
  array[2] int N_gp_ID_pred;
  array[Ngp] real gpx;
  array[2,N_gp_ID[1]] int ID_gp1;
  array[2,N_gp_ID[2]] int ID_gp2;
  array[2,N_gp_ID_pred[1]] int ID_pred_gp1;
  array[2,N_gp_ID_pred[2]] int ID_pred_gp2;
  matrix[N_pred, K] X_pred;
  matrix[N_pred, G] Z_pred;
}

parameters {
  matrix<lower= 0> [K, ncat] b;  // fixed effects regression coefficients
  matrix [G, ncat - 1] z;  // random effects regression coefficients
  real<lower= 0> sigma;
  real<lower= 0> sdgp;
  real<lower= 0> gpscale;
  array[2] matrix[Ngp, ncat-1] zgp;
}

transformed parameters {
  matrix[N, ncat] mu_logit;
  matrix[N, ncat] mu;
  matrix[K, ncat] beta;
  matrix[Ngp, ncat] gp_scale1;
  matrix[Ngp, ncat] gp_scale2;
  matrix[N, ncat] gp_mu;
  real lprior = 0;
  
  beta = log(b);

  // set up the Gaussian process
  gp_scale1[,1] = rep_vector(0.0, Ngp);
  gp_scale2[,1] = rep_vector(0.0, Ngp);
  
  {
    matrix[Ngp, Ngp] L;
    matrix[Ngp, Ngp] gpK = gp_exp_quad_cov(gpx, sdgp, gpscale + 1);

    // diagonal elements
    for (n in 1:Ngp) {
      gpK[n, n] = gpK[n, n] + 0.0001;
    }

    L = cholesky_decompose(gpK);
    gp_scale1[,2:4] = L * zgp[1];
    gp_scale2[,2:4] = L * zgp[2];
  }
  
  gp_mu[ID_gp1[1]] = gp_scale1[ID_gp1[2]];
  gp_mu[ID_gp2[1]] = gp_scale2[ID_gp2[2]];
  
  mu_logit = X * beta + gp_mu + sigma * Z * append_col(rep_vector(0.0, G),z);
  
  for(i in 1:N){
    mu[i] = to_row_vector(softmax(to_vector(mu_logit[i])));
  }
  for(i in 1: K){
    lprior += gamma_lpdf(b[i] | 1, 1);
  }
  
  lprior +=  std_normal_lpdf(sdgp)
    - 2 * normal_lccdf(0 | 0,1);
  lprior += std_normal_lpdf(gpscale);
  lprior += student_t_lpdf(sigma | 3, 0, 0.1) - student_t_lccdf(0 | 3,0,0.1);
}

model {
  target += lprior;
  for(i in 1:2){
    target += std_normal_lpdf(to_vector(zgp[i]));
  }
  for (i in 1:G){
    target += std_normal_lpdf(z[i]);
  }
  for(i in 1:N){
    target += categorical_lpmf(Y[i] | to_vector(mu[i]));
  }
}

generated quantities {
  matrix[N_pred,ncat] gp_mu_pred;
  matrix[N_pred,ncat] mu_pred;
  vector[N] log_lik;
  
  for(i in 1:N){
    log_lik[i] = categorical_lpmf(Y[i] | to_vector(mu[i]));
  }
  
  {
    gp_mu_pred[ID_pred_gp1[1]] = gp_scale1[ID_pred_gp1[2]];
    gp_mu_pred[ID_pred_gp2[1]] = gp_scale2[ID_pred_gp2[2]];

    matrix[N_pred, ncat] mu_pred_logit = X_pred * beta + gp_mu_pred+  Z_pred * append_col(rep_vector(0.0, G),z);
    
    for(i in 1:N_pred){
      mu_pred[i] = to_row_vector(softmax(to_vector(mu_pred_logit[i])));
    }
  }
}
