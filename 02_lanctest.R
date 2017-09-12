#Testrun of first 500 SNPs
#cut -d " " -f 1-500 output/theo_linear_c2_6.ps21.txt > ch21test.txt

##TBD - exclude subjects w/o glaucoma!
#Data has subject per row, so correct format

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

library(data.table)

for (chr in  c(c("_1","_2","_3","_4","_5","_6","_7"),c(8:22)))
{
#Load Local ancestry
dat <- fread(paste("output/theo_linear_c2",chr,".ps21.txt",sep=''), header=F,data.table=F)

#Load positions
positions <- fread(paste("output/theo_linear_c2",chr,".snpinfo.txt",sep=''), header=T,data.table=F)



#Make dataframes for each ancestry
dat_afr <- dat[africans,seq(2,dim(dat)[2],by=2)]
dat_eur <- dat[africans,seq(1,dim(dat)[2],by=2)]


#Get null proportion from ancestry
null_afr <-  anc3[africans,]$africa*2
null_eur <-  anc3[africans,]$europe*2
Site_of_Sample_Preparation <- anc3[africans,]$Site_of_Sample_Preparation

#Do regression. For the case only test, regression is the bias corrected (due to slope not being 1) version of the paired t test
    t2 <- lm(as.matrix(dat_afr) ~ null_afr + Site_of_Sample_Preparation )
    t2s <- summary(t2)
    extractor <- function(x){ return(x$coefficients[1,4]) }
#Get regression p-values
    regsums <- sapply(t2s,extractor)

#Results
    outdf <- data.frame(cbind(subset(positions,select=c(rs,pos,chr,maf)),regsums))
    write.table(outdf, paste("results/chr",chr,"_caseonly_samplesite.results",sep=''),row.names=F)

##Plot results

#Average ancestry per allele (Y axis)
    # avg_african <- apply(dat_afr,2,mean,na.rm=T)
    # avg_european <- apply(dat_eur,2,mean,na.rm=T)

    # mean_dif <- function(x,difval) {
        # mean(x - difval,na.rm=T)
    # }

    # pos_dif <- apply(dat_afr,2,mean_dif,null_afr)


#Note which values had significant p
    # sigvals <- which(regsums < .005)
    # sig_colors <- rep(afr_color,length(avg_african))
    # sig_colors[sigvals] <- "red"

#Plot
    # pdf(paste('chr',chr,'_onlyglaucoma.pdf',sep=''),7,7)

    # newpos <- round(positions$pos/10000,1)
    # plot(newpos,avg_african,type="l",lwd=3,col=afr_color,ylim=c(.1,1.8),xlab="Position (KB)", ylab="Average Ancestry")

    # points(newpos[sigvals],avg_african[sigvals],pch=16,col="red")
    # abline(h=mean(null_afr),lwd=2,col="black",lty=2)

    # lines(newpos,avg_european,type="l",lwd=3,col=eur_color)


    # plot(newpos,pos_dif ,type="l",lwd=3,col=afr_color,ylim=c(-.2,.2),xlab="Position (KB)", ylab="Average Ancestry Deviation")
    # points(newpos[sigvals],pos_dif[sigvals],pch=16,col="red")

    # dev.off()
}


###Extras



#What is the average sd per allele?
sd_african <- apply(dat_afr,2,sd,na.rm=T)
median(sd_african)
median(avg_african)



#lets look at one of these regions...
 head(positions[sigvals,])
 tail(positions[sigvals,])

positions[which(positions$pos >= 22000000 & positions$pos <= 22010000),]
regsums[22571] #Nominal assn (p=0.03) of rs1063192 on chr 9

regsums[51222] #Nominal assn (p=0.004) of rs735860 on chr 6



#Chromosomal level admixture proportions
admixes0 <- read.table('output/theo_linear_c210.admix.txt',header=F)
names(admixes0) <- c("eur","afr")
admixes2 <- data.frame(admixes0,subs)
admixes3 <- merge(admixes2,anc,by=c("FID","IID"),all.x=TRUE)
admixes <- admixes3[order(admixes3$order),]



t1 <- summary(lm(dat_afr[,1] ~ null_afr , data=dat_afr))
t.test(dat_afr[,1],dat_afr$null_afr,paired=T)

#Average ancestry per subject
avg_subj <- apply(dat_afr,1,mean)

mean(avg_subj)

cor.test(null_afr,avg_subj)
