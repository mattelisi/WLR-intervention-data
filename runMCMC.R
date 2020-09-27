rm(list=ls())
setwd('~/git_local/WLR-intervention-data/')

library(rstan)
options(mc.cores = parallel::detectCores()) # indicate stan to use multiple cores if available
rstan_options(auto_write = TRUE)

# happy sad

# full model
d_stan <- readRDS("./data_stan/biasHS_full.RDS")
m0 <- stan(file = "./stan_model/bGLMM_bias_WLR.stan", data = d_stan, iter = 2000, chains = 4)
saveRDS(m0, file="./stan_results/bGLMM_bias_WLR_results_full.RDS")

# full model with correction
d_stan <- readRDS("./data_stan/biasHS_full.RDS")
m0 <- stan(file = "./stan_model/bGLMM_bias_WLR_corrected.stan", data = d_stan, iter = 2000, chains = 4)
saveRDS(m0, file="./stan_results/bGLMM_bias_WLR_results_full_corrected.RDS")

# fear anger - full model
# full model with correction
d_stan <- readRDS("./data_stan/biasFA_full.RDS")
m0 <- stan(file = "./stan_model/bGLMM_bias_WLR_corrected_FA.stan", data = d_stan, iter = 2000, chains = 4)
saveRDS(m0, file="./stan_results/bGLMM_bias_WLR_results_full_corrected_FA.RDS")

