# prepare data for analysis of HS task
# clear workspace and set working directory
rm(list=ls())
setwd('~/git_local/WLR-intervention-data/')

dMerged <- read.table("./data_bias/bias_T123_merged.txt",sep=";",header=T)

# JM9 came to testing at T3 but we then discovered that had stopped attending the program
# I use their data only for baseline at T1
#dMerged <- dMerged[dMerged$id!="JM9",]
dMerged <- dMerged[-which(dMerged$id=="JM9" & dMerged$time!="T1"),]
d <- dMerged[dMerged$task=="HS",]



# ------------------------------------------------------------------------------------------------------------- #
### make also data for Stan model - HS task
d <- dMerged[dMerged$task=="HS",]

# use deviation (sum-to-zero) coding for face_gender
contrasts(d$face_gender) <- "contr.sum"
d$face_gender_dummy <- ifelse(d$face_gender=="F",1,-1)
d$id_n <- as.numeric(factor(d$id, labels=1:length(unique(d$id))))
d$time_n <- as.numeric(factor(d$time, labels=1:length(unique(d$time))))
tapply(d$time_n, d$time, mean) # sanity
d$group_n <- ifelse(d$group=="E",1,0)
tapply(d$group_n, d$group, mean) # sanity

unique(d$age)

unique(d$actor)
d$actor_dummy <- ifelse(d$actor==2,1,0)
tapply(d$actor_dummy, list(d$group,d$time,d$ethni),mean) # sanity

saveRDS(d,"./data_bias/HSdata.RDS")

#compile stan data list
d_stan <- list(N=nrow(d),
               J=length(unique(d$id)),
               id=d$id_n,
               rfg = d$refugee,
               face_gender = d$face_gender_dummy,
               actor = d$actor_dummy ,
               T2 = ifelse(d$time_n==2,1,0),
               T3 = ifelse(d$time_n==3,1,0),
               group = d$group_n,
               classification = d$classification,
               img_number = d$img_number,
               age = d$age
)
str(d_stan)

# write to file
saveRDS(d_stan, file="./data_stan/biasHS_full.RDS")
