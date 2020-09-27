data {
  int<lower=1> N;                         // number of observations
  int<lower=1> J;                         // number of subjects
  int<lower=0,upper=1> classification[N]; // response
  int<lower=1,upper=J> id[N];             // subject id
  vector[N] img_number;                   // img_number correspond to morph level
  int<lower=-1,upper=1> face_gender[N];   // contr.sum for face gender (F=1, M=-1)
  int<lower=0,upper=1> T2[N];             // time point 2
  int<lower=0,upper=1> T3[N];             // time point 3
  int<lower=0,upper=1> group[N];          // set to 1 for experimental (WLR)
  int<lower=0,upper=1> rfg[N];            // refugee status
  int<lower=0,upper=1> actor[N];          // flag == 1 for actor 2
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
  vector[2] beta_rfg;          // refugee
  vector[2] beta_group;        // group (C=0;E=1)
  
  matrix[2,2] beta_time;           // [locT2, locT3 refugees; locT2, locT3 controls]
  matrix[2,2] beta_time_groupE;    // [locT2, locT3 refugees; locT2, locT3 controls] experimental group
  
  matrix[2,2] beta2_time;           // [scaleT2, scaleT3 refugees; scaleT2, scaleT3 controls]
  matrix[2,2] beta2_time_groupE;    // [scaleT2, scaleT3 refugees; scaleT2, scaleT3 controls] experimental group
  
  // parameters for adjustment for Syrianc Control T3
  vector[2] beta_actor;
  vector[2] beta_actor_genderF;

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
  beta_rfg ~ normal(0, 1);
  to_vector(beta_time) ~ normal(0, 1);
  to_vector(beta2_time) ~ normal(0, 1);
  beta_group ~ normal(0, 1);
  to_vector(beta_time_groupE) ~ normal(0,1);
  to_vector(beta2_time_groupE) ~ normal(0,1);
  
  // prior for adjustment parameters
  // obtained for analysis of the other dataset
  
  // this is happy-sad
  beta_actor[1] ~ normal(-0.5687, 0.5022);
  beta_actor[2] ~ normal(0.3409, 0.6211);
  beta_actor_genderF[1] ~ normal(-0.0918, 0.6338);
  beta_actor_genderF[2] ~ normal(-0.3555, 0.6919);

  /*
  //fear - anger
  beta_actor[1] ~ normal(0.7434, 0.3111);
  beta_actor[2] ~ normal(0.6194, 0.4075);
  beta_actor_genderF[1] ~ normal(-0.4076, 0.4300);
  beta_actor_genderF[2] ~ normal(-0.4995, 0.5595);
  */

  lambda ~ beta(0.5, 15);    // lapse rate (~50% prior belief that lapse rate <=0.1)
  
  //likelihood
  for (i in 1:N){
    
    mu = beta[1] + u[1,id[i]] + face_gender[i]*beta_face_gender[1] + beta_rfg[1]*rfg[i] + group[i]*beta_group[1]
      + actor[i]*beta_actor[1] + actor[i]*((face_gender[i]+1)/2*beta_actor_genderF[1])
      + T2[i] * (beta_time[rfg[i]+1,1] + group[i]*beta_time_groupE[rfg[i]+1,1])
      + T3[i] * (beta_time[rfg[i]+1,2] + group[i]*beta_time_groupE[rfg[i]+1,2])
      + morphl_Z[i] * (beta[2] + u[2,id[i]] + face_gender[i]*beta_face_gender[2] + beta_rfg[2]*rfg[i]  + group[i]*beta_group[2]
        + actor[i]*beta_actor[2] + actor[i]*((face_gender[i]+1)/2*beta_actor_genderF[2])
        + T2[i] * (beta2_time[rfg[i]+1,1] + group[i]*beta2_time_groupE[rfg[i]+1,1])
        + T3[i] * (beta2_time[rfg[i]+1,2] + group[i]*beta2_time_groupE[rfg[i]+1,2]));
    
    target += bernoulli_lpmf(classification[i] | (1-lambda)*Phi_approx(mu)+lambda/2);
  }
}

generated quantities {
  vector[N] log_lik;
  real mu2; // conditional mean of the dependent variable
  for (i in 1:N) {
    mu2 = beta[1] + u[1,id[i]] + face_gender[i]*beta_face_gender[1] + beta_rfg[1]*rfg[i] + group[i]*beta_group[1]
      + actor[i]*beta_actor[1] + actor[i]*((face_gender[i]+1)/2*beta_actor_genderF[1])
      + T2[i] * (beta_time[rfg[i]+1,1] + group[i]*beta_time_groupE[rfg[i]+1,1])
      + T3[i] * (beta_time[rfg[i]+1,2] + group[i]*beta_time_groupE[rfg[i]+1,2])
      + morphl_Z[i] * (beta[2] + u[2,id[i]] + face_gender[i]*beta_face_gender[2] + beta_rfg[2]*rfg[i]  + group[i]*beta_group[2]
        + actor[i]*beta_actor[2] + actor[i]*((face_gender[i]+1)/2*beta_actor_genderF[2])
        + T2[i] * (beta2_time[rfg[i]+1,1] + group[i]*beta2_time_groupE[rfg[i]+1,1])
        + T3[i] * (beta2_time[rfg[i]+1,2] + group[i]*beta2_time_groupE[rfg[i]+1,2]));
        
    log_lik[i] = bernoulli_lpmf(classification[i] | (1-lambda)*Phi_approx(mu2)+lambda/2);
  }
}
