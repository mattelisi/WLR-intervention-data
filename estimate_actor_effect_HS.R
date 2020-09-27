# ------------------------------------------------------------------------------------------------------------- #
# This script run the additional analysis done to estiamate 
# the effect of the change of identity in the sitmulus for a subset of children at T3

rm(list=ls())
setwd('~/git_local/WLR-intervention-data/')
dMerged <- read.table("./data_bias/bias_T123_merged.txt",sep=";",header=T)
dMerged <- dMerged[dMerged$id!="JM9",]
d <- dMerged[dMerged$task=="HS",]

dag <- aggregate(refugee~id+group+ethni+time+face_gender+actor, d,mean)
aggregate(id~group+ethni+refugee+actor+face_gender, dag[dag$time=="T1",], length)
aggregate(id~group+ethni+refugee+actor+face_gender, dag[dag$time=="T2",], length)
aggregate(id~group+ethni+refugee+actor+face_gender, dag[dag$time=="T3",], length)
unique(d$id[d$actor==2])
# IM36 IM37 IM38 JM38 JM39 JM40 JM43 JM48 JM50
# these have actor 2 at T3
# these are all Syrian control, tested on 20th June 
# Note: these are from the first day of testing at the community center
# due to a mix-up in the scheduling they were mistakenly treated as new children 
# not part of the intervention

# set up extra databse
d_check <- read.table("./data_bias/additional_data/bias_all.txt",sep=";",header=T)
group_index <- d_check$time!= "T2" & (d_check$time!= "T3" | d_check$group== "N_T3")
d_check <- d_check[which(group_index),]
unique(d_check$time)
d_check <- d_check[which(d_check$task== "HS"),]
d_check <- d_check[which(d_check$ethni!= "B"),]
d_check <- d_check[which(d_check$ethni== "S"),]
d_check <- d_check[which(d_check$time!= "SSC"),]
length(unique(d_check$id))
d_check$refugee <- ifelse(d_check$ethni=="S",1,0)
dag_check <- aggregate(refugee~id+group+ethni+face_gender+actor+time+group, d_check, mean)
tab4 <- aggregate(id~actor+face_gender+ethni, dag_check, length)

tab4$ethni<- ifelse(tab4$ethni=="J","Jordanian","Syrian")
colnames(tab4)[1] <- "Actor identity"
colnames(tab4)[2] <- "Actor gender"
colnames(tab4)[3] <- "Nationality"
colnames(tab4)[4] <- "n. subjects"

with(dag_check, tapply(id, list(time,group),length))

# set deviation (sum-to-zero) coding for face_gender
d_check$face_gender_dummy <- ifelse(d_check$face_gender=="F",1,0)
d_check$id_n <- as.numeric(factor(d_check$id, labels=1:length(unique(d_check$id))))
d_check$group_n <- ifelse(d_check$group=="E",1,0)
tapply(d_check$group_n, d_check$group, mean) # sanity
d_check$actor_dummy <- ifelse(d_check$actor==2,1,0)

#compile stan data list
d_stan <- list(N=nrow(d_check),
               J=length(unique(d_check$id)),
               id=d_check$id_n,
               rfg = d_check$refugee,
               face_gender = d_check$face_gender_dummy,
               actor = d_check$actor_dummy,
               group = d_check$group_n,
               classification = d_check$classification,
               img_number = d_check$img_number 
)
str(d_stan)

library(rstan)
options(mc.cores = parallel::detectCores()) # indicate stan to use multiple cores if available
rstan_options(auto_write = TRUE)
#m0 <- stan(file = "./additional_data/bGLMM_bias_check_rfgOnly.stan", data = d_stan, iter = 3000, chains = 4)
#saveRDS(m0, file="./additional_data/bGLMM_bias_check_results_rfgOnly.RDS")

m0 <- readRDS("./data_bias/additional_data/bGLMM_bias_check_results_rfgOnly.RDS")

# check convergence
traceplot(m0, pars=c("beta","beta_face_gender","beta_actor","beta_actorXgender","lambda"))
#print(m0, pars=c("beta","beta_face_gender","beta_actor","beta_actorXgender","lambda","sigma_u"))

# extract psychometric parameters
m0@model_pars
beta1 <- extract(m0, pars=c("beta[1]"))$"beta[1]"
beta2 <- extract(m0, pars=c("beta[2]"))$"beta[2]"
#beta_rfg <- extract(m0, pars=c("beta_rfg"))$"beta_rfg"
beta_gender <- extract(m0, pars=c("beta_face_gender"))$"beta_face_gender"
beta_actor <- extract(m0, pars=c("beta_actor"))$"beta_actor"
beta_actorXgender <- extract(m0, pars=c("beta_actorXgender"))$"beta_actorXgender"

B0_F1 <- beta1 + beta_gender[,1] #+ beta_rfg[,1]
B0_F2 <- beta1 + beta_gender[,1] + beta_actor[,1] + beta_actorXgender[,1] #+ beta_rfg[,1]
B0_M1 <- beta1 # + beta_rfg[,1]
B0_M2 <- beta1 + beta_actor[,1] #+ beta_rfg[,1]

B1_F1 <- beta2 + beta_gender[,2] #+ beta_rfg[,2]
B1_F2 <- beta2 + beta_gender[,2] + beta_actor[,2] + beta_actorXgender[,2] #+ beta_rfg[,2]
B1_M1 <- beta2 # + beta_rfg[,2]
B1_M2 <- beta2 + beta_actor[,2] #+ beta_rfg[,2]

scl_val <- unique(abs(range((d_stan$img_number- 11) / 6.2)))

mu_F1 <- (-B0_F1 / B1_F1 )/scl_val 
mu_F2 <- (-B0_F2 / B1_F2 )/scl_val 
mu_M1 <- (-B0_M1 / B1_M1 )/scl_val 
mu_M2 <- (-B0_M2 / B1_M2 )/scl_val 

sig_F1 <- (1 / B1_F1 )/scl_val 
sig_F2 <- (1/ B1_F2 )/scl_val 
sig_M1 <- (1/ B1_M1 )/scl_val 
sig_M2 <- (1/ B1_M2 )/scl_val 

# posterior distribution and Gaussian approximation
#pdf("./additional_data/actor_eff_beta.pdf",width=6,height=6)
par(mfrow=c(2,2))
hist(beta_actor[,1], breaks=40, freq=F, col="grey",border="white",main="",xlab=expression(beta[3]))
curve(dnorm(x, mean = mean(beta_actor[,1]),sd=sd(beta_actor[,1])), col = 2,  lwd = 2, add = TRUE)

hist(beta_actor[,2], breaks=40, freq=F, col="grey",border="white", main="",xlab=expression(beta[6]))
curve(dnorm(x, mean = mean(beta_actor[,2]),sd=sd(beta_actor[,2])), col = 2,  lwd = 2, add = TRUE)

hist(beta_actorXgender[,1], breaks=40, freq=F, col="grey",border="white", main="",xlab=expression(beta[4]))
curve(dnorm(x, mean = mean(beta_actorXgender[,1]),sd=sd(beta_actorXgender[,1])), col = 2,  lwd = 2, add = TRUE)

hist(beta_actorXgender[,2], breaks=40, freq=F, col="grey",border="white", main="",xlab=expression(beta[7]))
curve(dnorm(x, mean = mean(beta_actorXgender[,2]),sd=sd(beta_actorXgender[,2])), col = 2,  lwd = 2, add = TRUE)
#dev.off()

# make dataframe to store informative prior for full model
prior_actor <- data.frame(parameter=c(rep("intercept",2),rep("slope",2)), 
                          stimulus = c("M2","F2","M2","F2"), 
                          Mean = c(mean(beta_actor[,1]), mean(beta_actorXgender[,1]), mean(beta_actor[,2]), mean(beta_actorXgender[,2])),
                          Std = c(sd(beta_actor[,1]), sd(beta_actorXgender[,1]), sd(beta_actor[,2]), sd(beta_actorXgender[,2])))
write.table(prior_actor, file="./data_bias/additional_data/prior_actor.txt",row.names=F,quote=F)

