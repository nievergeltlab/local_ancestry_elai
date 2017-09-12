library(coda)
library(stats)
library(data.table)

##Estimate teh number of switches for each subject on each chromosome


lat_color <- rgb(230,185,184,max=255)
eur_color <- rgb(147,205, 221,max=255)
afr_color <- rgb(179,162,199,max=255)
afr_color_sig <- rgb(255,162,199,max=255)


#Load overall admixture proportions (SNPweights AIMs panelk estimates)
anc <- read.table('ancestry_missing_set_to_ref2.predpc_oneweek.header',header=T)
phenotypes <- read.table("/mnt/sdb/genetics/ancestry_pipeline_elai/glaucoma_shiley_jun26.pheno",header=T,stringsAsFactors=F,na.strings=c("NA","#N/A"))

subs <- read.table('tf/theo2mil_v6B_5.fam',header=F)
names(subs)[c(1:2)] <- c("FID","IID")

subs$order <- 1:dim(subs)[1]

anc1 <- merge(anc,subs,by=c("FID","IID"),all.y=TRUE)
anc2 <- merge(phenotypes,anc1,by=c("FID"),all.y=TRUE)
anc3 <- anc2[order(anc2$order),]

#i'm going to look at glaucoma only subjects and add a covariate 
africans <- which(anc3$bestpop_oneweek=="aam" & anc3$Diagnosis == "Glaucoma") # anc3$Site_of_Sample_Preparation
#Need to load position info

efsizessum <- rep(0,length(africans))

for (chr in c(c(8:22))) #c("_1","_2","_3","_4","_5","_6","_7")
{

#Load Local ancestry
dat <- fread(paste("output/theo_linear_c2",chr,".ps21.txt",sep=''), header=F,data.table=F)

#Load positions
positions <- fread(paste("output/theo_linear_c2",chr,".snpinfo.txt",sep=''), header=T,data.table=F)

dat_afr <- dat[africans,seq(2,dim(dat)[2],by=2)]

#Make into a time series object
ts1 <- mcmc(t(dat_afr))

#Get effective swithc N for each
efsizes <- effectiveSize(ts1)

efsizessum <- efsizessum + efsizes
write.table(efsizes, paste('efsizes/',chr,"_.efsizes",sep=""))

}

mean(efsizessum) # Expected N swithces = 201.71~ 202