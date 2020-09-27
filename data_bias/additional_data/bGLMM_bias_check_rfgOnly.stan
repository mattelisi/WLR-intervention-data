data {
  int<lower=1> N;                         // number of observations
  int<lower=1> J;                         // number of subjects
  int<lower=0,upper=1> classification[N]; // response
  int<lower=1,upper=J> id[N];             // subject id
  vector[N] img_number;                   // img_number correspond to morph level
  int<lower=0,upper=1> face_gender[N];    // contr.sum for face gender
  int<lower=0,upper=1> actor[N];          //actor
}

transformed data {
  // transform predictors on unit scale
  vector[N] morphl_Z;
  morphl_Z = (img_number - 11) / 6.2;
}

parameters {
  // within subject parameters
  vector[2] beta;               // fixed-effects parameters
  vector<lower=0>[2] sigma_u;   // random effects standard deviations
  cholesky_factor_corr[2] L_u;  // L_u is the Choleski factor of the correlation matrix
  matrix[2,J] z_u;              // random effect matrix
  
  // group-level parameters    // all: [location, scale]
  vector[2] beta_face_gender;  // face
  vector[2] beta_actor;        // actor
  vector[2] beta_actorXgender; // interaction

  real<lower=0,upper=1> lambda; // lapse rate (bounded)
}

transformed parameters {
  matrix[2,J] u;
  u = diag_pre_multiply(sigma_u, L_u) * z_u; // use Cholesky to set correlation
}

model {
  real mu; // conditional mean of the dependent variable

  //priors
  L_u ~ lkj_corr_cholesky(2);   // LKJ prior for the correlation matrix
  to_vector(z_u) ~ normal(0,1); // before Cholesky, random effects are normal variates with SD=1
  sigma_u ~ cauchy(0, 1);       // SD of random effects 
  beta[1] ~ normal(0, 1);          // prior for fixed-effect intercept
  beta[2] ~ normal(2, 5);          // prior for fixed-effect slope (this is broader and positive)
  
  beta_face_gender ~ normal(0, 1);
  beta_actor ~ normal(0, 1);
  beta_actorXgender ~ normal(0, 1);

  lambda ~ beta(0.5, 15);    // lapse rate (~50% prior belief that lapse rate <=0.1)
  
  //likelihood
  for (i in 1:N){
    
    mu = beta[1] + u[1,id[i]] + face_gender[i]*beta_face_gender[1] + actor[i]*beta_actor[1] + face_gender[i]*actor[i]*beta_actorXgender[1] 
      + morphl_Z[i] * (beta[2] + u[2,id[i]] + face_gender[i]*beta_face_gender[2] + actor[i]*beta_actor[2] + face_gender[i]*actor[i]*beta_actorXgender[2]);
    
    target += bernoulli_lpmf(classification[i] | (1-lambda)*Phi_approx(mu)+lambda/2);
  }
}
