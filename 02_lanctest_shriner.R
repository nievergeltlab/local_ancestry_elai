
library(data.table)

#Plot colors
 eur_color <- rgb(147,205, 221,max=255)
 afr_color <- rgb(179,162,199,max=255)
 afr_color_sig <- rgb(255,162,199,max=255)

#Load overall admixture proportions (SNPweights AIMs panel estimates) for case data
 anc <- read.table('/mnt/sdb/genetics/ancestry_pipeline_elai/ancestry_missing_set_to_ref2.predpc_oneweek.header',header=T)
 phenotypes <- read.table("/mnt/sdb/genetics/ancestry_pipeline_elai/glaucoma_shiley_jun26.pheno",header=T,stringsAsFactors=F,na.strings=c("NA","#N/A"))

 subs <- read.table('/mnt/sdb/genetics/ancestry_pipeline_elai/tf/theo2mil_v6B_5.fam',header=F)
 names(subs)[c(1:2)] <- c("FID","IID")

 subs$order <- 1:dim(subs)[1]

 anc1 <- merge(anc,subs,by=c("FID","IID"),all.y=TRUE)
 anc2 <- merge(phenotypes,anc1,by=c("FID"),all.y=TRUE)
 anc3 <- anc2[order(anc2$order),]

#Subset to glaucoma only, black subjects
 africans <- which(anc3$bestpop_oneweek=="aam" & anc3$Diagnosis == "Glaucoma")

#get a covariate sheet
 anc3_cov <- anc3[africans,]

#Load overall admixture proportions (SNPweights AIMs panel estimates) for CONTROL data
 anc <- read.table('/mnt/sdb/genetics/gradyanc/GTPC/ancestry/Adam_final_GTP_pheno.predpc_oneweek.header',header=T)

 subs <- read.table('gtp_famsplit10',header=F)
 names(subs)[c(1:2)] <- c("FID","IID")

 subs$order <- 1:dim(subs)[1]

 anc2 <- merge(anc,subs,by=c("IID"),all.y=TRUE)
 anc3_cont <- anc2[order(anc2$order),]

#Subset to black subjects
 africans_cont <- which(anc3_cont$bestpop_oneweek=="aam" ) 

 #get a covariate sheet
 anc3_cov_cont <- anc3_cont[africans_cont ,]

#Combine case/control covariates sheets (PCs) for analysis
 ancov <- rbind(anc3_cov[,c("PC1","PC2","PC3","PC4","PC5")], anc3_cov_cont[,c("PC1","PC2","PC3","PC4","PC5")])


for (chr in  c(c("_1","_2","_3","_4","_5","_6","_7"),c(8:22)))
{
##Load positions for glaucoma data
 positions <- fread(paste("/mnt/sdb/genetics/ancestry_pipeline_elai/output/theo_linear_c2",chr,".snpinfo.txt",sep=''), header=T,data.table=F)

#Load Local ancestry for Glaucoma data
 dat <- fread(paste("/mnt/sdb/genetics/ancestry_pipeline_elai/output/theo_linear_c2",chr,".ps21.txt",sep=''), header=F,data.table=F)

#Make dataframes for each ancestry class
 dat_afr <- dat[africans,seq(2,dim(dat)[2],by=2)]
 dat_eur <- dat[africans,seq(1,dim(dat)[2],by=2)]

##Load positions for control data
 pos_cont <- fread(paste("output/gtp_realquad_07",chr,".snpinfo.txt",sep=''), header=T,data.table=F)

#Load Local ancestry for control data
 dat_cont <- fread(paste("output/gtp_realquad_07",chr,".ps21.txt",sep=''), header=F,data.table=F)

#Make dataframes for each ancestry
 dat_cont_afr <- dat_cont[africans_cont,seq(2,dim(dat_cont)[2],by=2)]
 dat_cont_eur <- dat_cont[africans_cont,seq(1,dim(dat_cont)[2],by=2)]


##Identify markers where calls were made for both cases/controls
 pos_overlap <- positions[which(positions$rs %in% pos_cont$rs),]
 pos_cont_overlap <- pos_cont[which(pos_cont$rs %in% pos_overlap$rs),]
 
 #check that this is correct
 print(table(pos_overlap$rs == pos_cont_overlap$rs))

#Subset case and controls to only theseoverlapping positions
 dat_afr2 <- dat_afr[,which(positions$rs %in% pos_cont$rs)]

 dat_cont_afr2 <- dat_cont_afr[,which(pos_cont$rs %in% pos_overlap$rs)]

#Reset column names so that data can stack
 names(dat_afr2) <- c(1:dim(dat_afr2)[2])
 names(dat_cont_afr2) <- c(1:dim(dat_cont_afr2)[2])

#Stack case/control data frames
 dat_afr3 <- rbind(dat_afr2,dat_cont_afr2)

 
#Do regression of  case status on local ancestry 
 extractor <- function(x){ return(x$coefficients[2,4]) }
 glmf <- function(x)
 {
  t2 <- glm(c(rep(1,dim(dat_afr2)[1]),rep(0,dim(dat_cont_afr2)[1])) ~ x + ancov$PC1 + ancov$PC2 + ancov$PC3,family='binomial' )
  t2s <- summary(t2)
  return(extractor(t2s))
 }
 
#Get regression p-values
    regsums <- apply(as.matrix(dat_afr3),2,glmf)

#Results
    outdf <- data.frame(cbind(subset(pos_overlap,select=c(rs,pos,chr,maf)),regsums))
    write.table(outdf, paste("results/chr",chr,"_caseonly_samplesite.results",sep=''),row.names=F)

##Plot results

#Average ancestry per allele (Y axis)
     avg_african_case <- apply(dat_afr2,2,mean,na.rm=T)
     avg_african_cont <- apply(dat_cont_afr2,2,mean,na.rm=T)
     
    # avg_european <- apply(dat_eur,2,mean,na.rm=T)

    # mean_dif <- function(x,difval) {
        # mean(x - difval,na.rm=T)
    # }

    # pos_dif <- apply(dat_afr,2,mean_dif,null_afr)
     pos_dif <- avg_african_case - avg_african_cont

#Note which values had significant p
     sigvals <- which(regsums < .005)
     sig_colors <- rep(afr_color,length(avg_african_case ))
    #sig_colors[sigvals] <- "red"

#Plot
     pdf(paste('results/chr',chr,'_withcontglaucoma.pdf',sep=''),7,7)

     newpos <- round(pos_overlap$pos/10000,1)
     plot(newpos,avg_african_case,type="l",lwd=3,col=afr_color,ylim=c(.1,1.8),xlab="Position (KB)", ylab="Average Ancestry")
     lines(newpos,avg_african_cont,type="l",lwd=3,col=eur_color)

      
     points(newpos[sigvals],avg_african[sigvals],pch=16,col="red")
    # abline(h=mean(null_afr),lwd=2,col="black",lty=2)

     #lines(newpos,avg_european,type="l",lwd=3,col=eur_color)


     plot(newpos,pos_dif ,type="l",lwd=3,col=afr_color,ylim=c(-.2,.2),xlab="Position (KB)", ylab="Average Ancestry Deviation")
     points(newpos[sigvals],pos_dif[sigvals],pch=16,col="red")

     dev.off()
     inc=inc+1
}


###Extras





#What is the average sd per allele?
sd_african <- apply(dat_afr,2,sd,na.rm=T)
median(sd_african)
median(avg_african)


#View top loci:
head(outdf[order(outdf$regsums),]

#lets look at one of these regions...
 head(pos_overlap[sigvals,])
 tail(positions[sigvals,])

pos_overlap[which(pos_overlap$pos >= 22000000 & positions$pos <= 22010000),]
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
