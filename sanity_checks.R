setwd('~/git_local/WLR-intervention-data/')
dMerged <- read.table("./data_bias/bias_T123_merged.txt",sep=";",header=T)

dMerged <- dMerged[dMerged$id!="JM9",]
d <- dMerged[dMerged$task=="HS",]

d2 <- dMerged[dMerged$id=="JM9",]


# Hi Matteo,
# Yes, that’s what I’ve got as well. She wasn’t tested at t2 because her mother stopped reading and then she came to the centre at t3 and we realized it’s her while she was doing one of the tasks for the new cohort. I don’t have psychopathology data on her for t2 or t3.
# 
# 
# From: Matteo Lisi <m.lisi@qmul.ac.uk>
#   Date: Monday, November 11, 2019 at 12:37 PM
# To: Julia Michalek <j.michalek@qmul.ac.uk>
#   Subject: Re: bias + questionnaire results
# 
# Hi Julia
# 
# I am doing the analyses for the bias across timepoints. There is one subject "JM9" that according to the data I have was tested at T1 and at T3, but not at T2; when you have a moment could you check whether that's the case or I am missing something?
