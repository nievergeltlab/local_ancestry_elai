library(coda)
library(stats)
library(data.table)

##Estimate teh number of switches for each subject on each chromosome

study="shiley_350"

lat_color <- rgb(230,185,184,max=255)
eur_color <- rgb(147,205, 221,max=255)
afr_color <- rgb(179,162,199,max=255)
afr_color_sig <- rgb(255,162,199,max=255)

study="adages_cleaner"

#Load overall admixture proportions (SNPweights AIMs panelk estimates)
anc <- read.table('ancestry/adages_cleaner.predpc_oneweek.header',header=T)
ancx <- read.table('results/adages_cleaner_global.ancestry',header=T)

anct <- merge(anc,ancx,by=c("FID","IID"))
summary(lm(africa ~ global_african,data=anct)) #High correspondence between veriones



phenotypes <- read.csv("jerry_glaucoma.csv",header=T,stringsAsFactors=F,na.strings=c("NA","#N/A"))
names(phenotypes)[1] <- "FID"

subs <- read.table('adages_cleaner.fam',header=F)
names(subs)[c(1:2)] <- c("FID","IID")

subs$order <- 1:dim(subs)[1]

anc1 <- merge(anc,subs,by=c("FID","IID"),all.y=TRUE)
anc2 <- merge(phenotypes,anc1,by=c("FID"),all.y=TRUE)
anc3 <- anc2[order(anc2$order),]

#i'm going to look at glaucoma only subjects and add a covariate 
africans <- which(anc3$bestpop_oneweek=="aam" & anc3$expected_group == "GLAUCOMA") # anc3$Site_of_Sample_Preparation
#africans <- which(anc3$best >= 0.15 & anc3$expected_group == "GLAUCOMA") # anc3$Site_of_Sample_Preparation

afrcov <- anc3[africans,]

#Need to load position info


efsizessum <- rep(0,length(africans))

for (chr in c(c(1:22))) 
{
print(chr)
#Load Local ancestry
 load(paste("output/",study,"_v4e_",chr,".ps21.header.aam.R",sep=""))
 dat <- data_combined

#Load positions
 positions=paste("output/",study,"_v4e_00_", chr,".snpinfo.txt",sep="") 
 posa <- fread(positions,data.table=F)
 pos_subset <- subset(posa,select=c(rs,pos,chr,maf))
 

#Make dataframes for each ancestry
 dat_afr <- dat[africans,-c(1:2)] #[africans,seq(2,dim(dat)[2],by=2)]
 dat_afr <- dat_afr[, which(names(dat_afr) %in% pos_subset$rs)]
 dat_afr <- round(dat_afr,0)
 
#Make into a time series object
 ts1 <- mcmc(t(dat_afr))

#Impute missing values
ts2 <- t(apply(ts1,1,function(x) {
  if(is.numeric(x)) ifelse(is.na(x),median(x,na.rm=T),x) else x}))
 

#Get effective swithc N for each
efsizes <- effectiveSize(ts2)

efsizessum <- efsizessum + efsizes
write.table(efsizes, paste('results/',chr,"_.efsizes",sep=""))

}

mean(efsizessum) # Expected N swithces = 201.71~ 202
