args <- commandArgs(trailingOnly = TRUE)
study <- args[1]
chr <- args[2] 
#anc_file <- args[3]
#pheno_file <- args[4]


##TBD - exclude subjects w/o glaucoma!
#Data has subject per row, so correct format

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

library(data.table)

for (chr in  c(c(1:22)))
{
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
 #dat_eur <- dat[africans,seq(1,dim(dat)[2],by=2)]



#Get null proportion from ancestry
 null_afr <-  anc3[africans,]$africa*2
 #null_afr <-  anc3[africans,]$global_african*2
 
 #null_eur <-  anc3[africans,]$europe*2
 #Site_of_Sample_Preparation <- anc3[africans,]$Site_of_Sample_Preparation

#Do regression. For the case only test, regression is the bias corrected (due to slope not being 1) version of the paired t test. The intercept is the mean difference needed.
#Cutting out subject IDs from dimension of data


 
 t2 <- lm(as.matrix(dat_afr) ~ null_afr + afrcov$PC1 + afrcov$PC2 + afrcov$PC3 + afrcov$PC4 + afrcov$PC5 )
 t2s <- summary(t2)
 extractor <- function(x){ return(x$coefficients[1,4]) }
#Get regression p-values
 regsums <- data.frame(sapply(t2s,extractor))
 names(regsums)[1] <- "p"
 regsums$rs <- gsub("Response ", "", row.names(regsums))
 
#Results

 outdf <- merge(pos_subset,regsums,by="rs")

 write.table(outdf[order(outdf$p),], paste("results/chr",chr,"_caseonly.results",sep=''),row.names=F)

##Plot results

#Average ancestry per allele (Y axis)
 avg_african <- apply(dat_afr,2,mean,na.rm=T)


 #avg_european <- apply(dat_eur,2,mean,na.rm=T)

 mean_dif <- function(x,difval) {
     mean(x - difval,na.rm=T)
 }

 #The null proportion is useless because it doesnt include all markers!
 pos_dif <- apply(dat_afr,2,mean_dif,null_afr)


#Note which values had significant p
 sigvals <- which(regsums < .00025)
 sig_colors <- rep(afr_color,length(avg_african))
 sig_colors[sigvals] <- "red"

#Plot
 pdf(paste('results/chr',chr,'_caseonly.pdf',sep=''),7,7)

 newpos <- round(pos_subset$pos/10000,1)
 plot(newpos,avg_african,type="l",lwd=3,col=afr_color,ylim=c(.1,1.8),xlab="Position (KB)", ylab="Average Ancestry")

 points(newpos[sigvals],avg_african[sigvals],pch=16,col="red")
 abline(h=mean(null_afr),lwd=2,col="black",lty=2)

 #lines(newpos,avg_european,type="l",lwd=3,col=eur_color)


 plot(newpos,pos_dif ,type="l",lwd=3,col=afr_color,ylim=c(-.2,.2),xlab="Position (KB)", ylab="Average Ancestry Deviation")
 points(newpos[sigvals],pos_dif[sigvals],pch=16,col="red")

 dev.off()
}

