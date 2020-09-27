# make a plot of psychometric functions
# with also the data

setwd('~/git_local/WLR-intervention-data/')
library(ggplot2)
source("/home/matteo/sync/miscR/miscFunctions.R")
library(rstan)
require(rethinking) 

# -------------------------------------------------------------------------------------------------- #
# load data
dMerged <- read.table("./data_bias/bias_T123_merged.txt",sep=";",header=T)
dMerged <- dMerged[dMerged$id!="JM9",]
d <- dMerged[dMerged$task=="FA",]

# -------------------------------------------------------------------------------------------------- #
# load fitted model
m0 <- readRDS("./stan_results/bGLMM_bias_WLR_results_full_corrected_FA.RDS")

# -------------------------------------------------------------------------------------------------- #
# plot of raw data
d$x_bin <- cut(d$img_number,6)
dag <- aggregate(cbind(classification,img_number) ~ group+x_bin+ethni+time + id, d, mean)
dag2 <- aggregate(cbind(classification,img_number) ~ group+x_bin+ethni+time, dag, mean)
dag2$se <- aggregate(classification ~ group+x_bin+ethni+time, dag, bootMeanSE)$classification

x_breaks <- sort(unique(d$img_number))
x_labels <- c("100%\nanger","","","","","","","","","","neutral","","","","","","","","","","100%\nfear")

ggplot(dag2, aes(x=img_number,y=classification, color=group))+geom_vline(xintercept=11,color="grey",size=0.2,lty=2)+geom_hline(yintercept=0.5,color="grey",size=0.2,lty=2)+geom_errorbar(data=dag2,aes(ymin=classification-se,ymax=classification+se),width=0)+geom_point(data=dag2,pch=21)+nice_theme+scale_x_continuous(breaks= x_breaks, labels=x_labels,limits=range(d$img_number))+labs(x="morph level",y="choice probability")+scale_color_manual(values=c("dark grey","black"),name=" ")+facet_grid(ethni~time)+geom_line()+ theme(panel.spacing.x = unit(1, "lines"))

# -------------------------------------------------------------------------------------------------- #
# compute psychometric functions
beta1 <- extract(m0, pars=c("beta[1]"))$"beta[1]"
beta2 <- extract(m0, pars=c("beta[2]"))$"beta[2]"
beta_group <- extract(m0, pars=c("beta_group"))$"beta_group"
beta_time <- extract(m0, pars=c("beta_time"))$"beta_time"
beta_time_groupE <- extract(m0, pars=c("beta_time_groupE"))$"beta_time_groupE"
beta2_time <- extract(m0, pars=c("beta2_time"))$"beta2_time"
beta2_time_groupE <- extract(m0, pars=c("beta2_time_groupE"))$"beta2_time_groupE"
beta_rfg <- extract(m0, pars="beta_rfg")$"beta_rfg"
lambda <- extract(m0, pars="lambda")$"lambda"

# location
B0_J_C_T1 <- beta1
B0_J_C_T2 <- beta1 + beta_time[,1,1]
B0_J_C_T3 <- beta1 + beta_time[,1,2]
B0_J_E_T1 <- beta1 + beta_group[,1]
B0_J_E_T2 <- beta1 + beta_group[,1] + beta_time[,1,1] + beta_time_groupE[,1,1]
B0_J_E_T3 <- beta1 + beta_group[,1] + beta_time[,1,2] + beta_time_groupE[,1,2]
B0_S_C_T1 <- beta1 + beta_rfg[,1]
B0_S_C_T2 <- beta1 + beta_rfg[,1] + beta_time[,2,1]
B0_S_C_T3 <- beta1 + beta_rfg[,1] + beta_time[,2,2]
B0_S_E_T1 <- beta1 + beta_rfg[,1] + beta_group[,1]
B0_S_E_T2 <- beta1 + beta_rfg[,1] + beta_group[,1] + beta_time[,2,1] + beta_time_groupE[,2,1]
B0_S_E_T3 <- beta1 + beta_rfg[,1] + beta_group[,1] + beta_time[,2,2] + beta_time_groupE[,2,2]

# scale
B1_J_C_T1 <- beta2
B1_J_C_T2 <- beta2 + beta2_time[,1,1]
B1_J_C_T3 <- beta2 + beta2_time[,1,2]
B1_J_E_T1 <- beta2 + beta_group[,2]
B1_J_E_T2 <- beta2 + beta_group[,2] + beta2_time[,1,1] + beta2_time_groupE[,1,1]
B1_J_E_T3 <- beta2 + beta_group[,2] + beta2_time[,1,2] + beta2_time_groupE[,1,2]
B1_S_C_T1 <- beta2 + beta_rfg[,2]
B1_S_C_T2 <- beta2 + beta_rfg[,2] + beta2_time[,2,1]
B1_S_C_T3 <- beta2 + beta_rfg[,2] + beta2_time[,2,2]
B1_S_E_T1 <- beta2 + beta_rfg[,2] + beta_group[,2]
B1_S_E_T2 <- beta2 + beta_rfg[,2] + beta_group[,2]  + beta2_time[,2,1] + beta2_time_groupE[,2,1]
B1_S_E_T3 <- beta2 + beta_rfg[,2] + beta_group[,2]  + beta2_time[,2,2] + beta2_time_groupE[,2,2]

# # transform in prop. of morphing from neutral
# #scl_val <- unique(abs(range((d_stan$img_number- 11) / 6.2)))
# 
# # mu
# mu_J_C_T1 <- (-B0_J_C_T1 / B1_J_C_T1 )#/scl_val 
# mu_J_C_T2 <- (-B0_J_C_T2 / B1_J_C_T2 )#/scl_val 
# mu_J_C_T3 <- (-B0_J_C_T3 / B1_J_C_T3 )#/scl_val 
# 
# mu_J_E_T1 <- (-B0_J_E_T1 / B1_J_E_T1 )#/scl_val 
# mu_J_E_T2 <- (-B0_J_E_T2 / B1_J_E_T2 )#/scl_val 
# mu_J_E_T3 <- (-B0_J_E_T3 / B1_J_E_T3 )#/scl_val 
# 
# mu_S_C_T1 <- (-B0_S_C_T1 / B1_S_C_T1 )#/scl_val 
# mu_S_C_T2 <- (-B0_S_C_T2 / B1_S_C_T2 )#/scl_val 
# mu_S_C_T3 <- (-B0_S_C_T3 / B1_S_C_T3 )#/scl_val 
# 
# mu_S_E_T1 <- (-B0_S_E_T1 / B1_S_E_T1 )#/scl_val 
# mu_S_E_T2 <- (-B0_S_E_T2 / B1_S_E_T2 )#/scl_val 
# mu_S_E_T3 <- (-B0_S_E_T3 / B1_S_E_T3 )#/scl_val 
# 
# # sigma
# sigma_J_C_T1 <- (1 / B1_J_C_T1 )#/scl_val  
# sigma_J_C_T2 <- (1 / B1_J_C_T2 )#/scl_val 
# sigma_J_C_T3 <- (1 / B1_J_C_T3 )#/scl_val 
# 
# sigma_J_E_T1 <- (1 / B1_J_E_T1 )#/scl_val 
# sigma_J_E_T2 <- (1 / B1_J_E_T2 )#/scl_val 
# sigma_J_E_T3 <- (1 / B1_J_E_T3 )#/scl_val 
# 
# sigma_S_C_T1 <- (1 / B1_S_C_T1 )#/scl_val 
# sigma_S_C_T2 <- (1 / B1_S_C_T2 )#/scl_val 
# sigma_S_C_T3 <- (1 / B1_S_C_T3 )#/scl_val 
# 
# sigma_S_E_T1 <- (1 / B1_S_E_T1 )#/scl_val 
# sigma_S_E_T2 <- (1 / B1_S_E_T2 )#/scl_val 
# sigma_S_E_T3 <- (1 / B1_S_E_T3 )#/scl_val 

# -------------------------------------------------------------------------------------------------- #
# build dataframe with predictions
x_imgn <- seq(1,21,length.out=50)
x <- (x_imgn- 11) / 6.2
alpha <- 0.05

# OK!!!!!
# eval(parse(text="hist(B0_J_C_T1)"))

ethni_i <- c("J","S")
group_j <- c("C","E")
time_t <- c("T1","T2","T3")
d_psy <- {}

for(i in ethni_i){
  for(j in group_j){
    for(t in time_t){
      B0 <- eval(parse(text=paste("B0",i,j,t,sep="_")))
      B1 <- eval(parse(text=paste("B1",i,j,t,sep="_")))
      P <- (1-t(repmat(t(lambda),50,1))) * pnorm(t(repmat(t(B0),50,1))+outer(B1,x)) + t(repmat(t(lambda),50,1))/2
      dP <- data.frame(x=x, img_number = x_imgn, P=apply(P,2,mean),
                       P_ub = apply(P,2, function(x) quantile(x, probs = 1 - alpha/2)),
                       P_lb = apply(P,2, function(x) quantile(x, probs = alpha/2)),
                       group = j,
                       ethni = i,
                       time = t)
      d_psy <- rbind(d_psy, dP)
    }
  }
}

# -------------------------------------------------------------------------------------------------- #
# plot together with raw data

d$x_bin <- cut(d$img_number,6)
dag <- aggregate(cbind(classification,img_number) ~ group+x_bin+ethni+time + id, d, mean)
dag2 <- aggregate(cbind(classification,img_number) ~ group+x_bin+ethni+time, dag, mean)
dag2$se <- aggregate(classification ~ group+x_bin+ethni+time, dag, bootMeanSE)$classification

x_breaks <- sort(unique(d$img_number))
x_labels <- c("100%\nanger","","","","","","","","","","neutral","","","","","","","","","","100%\nfear")
d_psy$classification <- 0

levels(dag2$group) <- c("control","WLR")
levels(dag2$ethni) <- c("Jordanian","Syrian")
levels(d_psy$group) <- c("control","WLR")
levels(d_psy$ethni) <-c("Jordanian","Syrian")

plotFA <- ggplot(dag2[dag2$time!="T3",], aes(x=img_number,y=classification, color=group, fill=group))+geom_vline(xintercept=11,color="grey",size=0.2,lty=2)+geom_hline(yintercept=0.5,color="grey",size=0.2,lty=2)+geom_errorbar(data=dag2[dag2$time!="T3",],aes(ymin=classification-se,ymax=classification+se),width=0)+geom_point(data=dag2[dag2$time!="T3",],pch=21,fill=NA)+nice_theme+scale_x_continuous(breaks= x_breaks, labels=x_labels,limits=range(d$img_number))+labs(x="morph level",y="choice probability")+scale_color_manual(values=c("dark grey","dark green"),name=" ")+scale_fill_manual(values=c("dark grey","dark green"),name=" ")+facet_grid(ethni~time)+ theme(panel.spacing.x = unit(1, "lines"))+geom_line(data=d_psy[d_psy$time!="T3",],aes(y=P),size=0.8)+geom_ribbon(data=d_psy[d_psy$time!="T3",],aes(ymin=P_lb, ymax=P_ub), alpha = 0.3, size = 0, color=NA)#+ggtitle("Fear-Anger task")

ggsave("./figures/psy_FA.pdf",plotFA, width=5,height=2.5)


# -------------------------------------------------------------------------------------------------- #
# transform in prop. of morphing from neutral
scl_val <- unique(abs(range((d$img_number- 11) / 6.2)))

# mu
mu_J_C_T1 <- (-B0_J_C_T1 / B1_J_C_T1 )/scl_val
mu_J_C_T2 <- (-B0_J_C_T2 / B1_J_C_T2 )/scl_val
mu_J_C_T3 <- (-B0_J_C_T3 / B1_J_C_T3 )/scl_val

mu_J_E_T1 <- (-B0_J_E_T1 / B1_J_E_T1 )/scl_val
mu_J_E_T2 <- (-B0_J_E_T2 / B1_J_E_T2 )/scl_val
mu_J_E_T3 <- (-B0_J_E_T3 / B1_J_E_T3 )/scl_val

mu_S_C_T1 <- (-B0_S_C_T1 / B1_S_C_T1 )/scl_val
mu_S_C_T2 <- (-B0_S_C_T2 / B1_S_C_T2 )/scl_val
mu_S_C_T3 <- (-B0_S_C_T3 / B1_S_C_T3 )/scl_val

mu_S_E_T1 <- (-B0_S_E_T1 / B1_S_E_T1 )/scl_val
mu_S_E_T2 <- (-B0_S_E_T2 / B1_S_E_T2 )/scl_val
mu_S_E_T3 <- (-B0_S_E_T3 / B1_S_E_T3 )/scl_val

d_S_C <- mu_S_C_T2 - mu_S_C_T1  
d_S_E <- mu_S_E_T2 - mu_S_E_T1 
d_J_C <- mu_J_C_T2 - mu_J_C_T1 
d_J_E <- mu_J_E_T2 - mu_J_E_T1 

# lb_hpdi <- function(x) HPDI(unname(x), prob=0.95)[1]
# ub_hpdi <- function(x) HPDI(unname(x), prob=0.95)[2]
lb_hpdi <- function(x) quantile(x, probs = alpha/2)
ub_hpdi <- function(x) quantile(x, probs = 1 - alpha/2)
d_diff <- data.frame(ethni=c("Syrian","Syrian","Jordanian","Jordanian"),
                     group =  c("control","WLR", "control","WLR"),
                     delta = c(mean(d_S_C), mean(d_S_E), mean(d_J_C), mean(d_J_E)),
                     lb = c(lb_hpdi(d_S_C), lb_hpdi(d_S_E), lb_hpdi(d_J_C), lb_hpdi(d_J_E)),
                     ub = c(ub_hpdi(d_S_C), ub_hpdi(d_S_E), ub_hpdi(d_J_C), ub_hpdi(d_J_E))
)

diff_plot <- ggplot(d_diff, aes(x=delta*100, xmin=lb*100, xmax=ub*100, y=group, color=group))+facet_grid(ethni~.)+geom_vline(xintercept=0,lty=2,size=0.5)+geom_point(size=3)+geom_errorbarh(height=0.2)+nice_theme+scale_color_manual(values=c("dark grey","dark green"),name=" ",guide=F)+labs(y="", x="PSE difference from T1\n[% morphing]")+coord_cartesian(xlim=c(-25,25))
ggsave("./figures/CIdiff_FA.pdf",diff_plot, width=2.5,height=2)

