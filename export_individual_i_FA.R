# anger-fear WLR program
# save individual estimates

# clear workspace and set working directory
rm(list=ls())
setwd('~/git_local/WLR-intervention-data/')

# stan
library(rstan)
#options(mc.cores = parallel::detectCores()) # indicate stan to use multiple cores if available
#rstan_options(auto_write = TRUE)

# load data
d <- readRDS("./data_bias/FAdata.RDS")

# # run sampling (takes a while..)
#m0 <- stan(file = "./stan_model/bGLMM_bias_WLR.stan", data = d_stan, iter = 4000, chains = 4)
#saveRDS(m0, file="./stan_results/bGLMM_bias_WLR_results1.RDS")

# load fitted model
m0 <- readRDS("./stan_results/bGLMM_bias_WLR_results_full_corrected_FA.RDS")
#pairs(m0, pars=c("beta_time","beta_time_groupE"))
m0@model_pars
print(m0, pars=c("beta","beta_rfg","beta_group","beta_time","beta_time_groupE","beta2_time","beta2_time_groupE","beta_face_gender","lambda"),digits=2,probs=c(0.05/2, 1-0.05/2))

# ------------------------------------------------------------------------------------------------------------- #
# compute group-level parameters
beta1 <- extract(m0, pars=c("beta[1]"))$"beta[1]"
beta2 <- extract(m0, pars=c("beta[2]"))$"beta[2]"
beta_group <- extract(m0, pars=c("beta_group"))$"beta_group"
beta_time <- extract(m0, pars=c("beta_time"))$"beta_time"
beta_time_groupE <- extract(m0, pars=c("beta_time_groupE"))$"beta_time_groupE"
beta2_time <- extract(m0, pars=c("beta2_time"))$"beta2_time"
beta2_time_groupE <- extract(m0, pars=c("beta2_time_groupE"))$"beta2_time_groupE"
beta_rfg <- extract(m0, pars="beta_rfg")$"beta_rfg"
lambda <- extract(m0, pars="lambda")$"lambda"

# individual random effects
U <- extract(m0, pars=c("u"))$"u"

# ------------------------------------------------------------------------------------------------------------- #
str(d)
d_id <- aggregate(id_n ~ id + refugee + ethni + group  + group_n, d, mean)

d_id$mu <- NA
d_id$sigma <- NA
for(i in 1:nrow(d_id)){
  n_id <- d_id$id_n[i]
  B_0 <- beta1 + d_id$refugee[i]*beta_rfg[,1] + d_id$group_n[i]*beta_group[,1] + U[,1,n_id]
  B_1 <- beta2 + d_id$refugee[i]*beta_rfg[,2] + d_id$group_n[i]*beta_group[,2] + U[,2,n_id]
  d_id$mu[i] <- -mean(B_0)/mean(B_1)
  d_id$sigma[i] <- 1/mean(B_1)
}

plot(d_id$mu, d_id$sigma)
str(d_id)

write.table(d_id, file="./individualT1_FA.txt",sep=";",quote=F,row.names=F)

# ------------------------------------------------------------------------------------------------------------- #

