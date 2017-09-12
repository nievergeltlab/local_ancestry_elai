
library(data.table)

#Plot colors
 eur_color <- rgb(147,205, 221,max=255)
 afr_color <- rgb(179,162,199,max=255)
 afr_color_sig <- rgb(255,162,199,max=255)

#Load overall admixture proportions (SNPweights AIMs panel estimates) for case data
 anc <- read.table('/mnt/sdb/genetics/ancestry_pipeline_elai/ancestry_missing_set_to_ref2.predpc_oneweek.header',header=T,stringsAsFactors=F)
 
 
 phenotypes <- read.table("/mnt/sdb/genetics/ancestry_pipeline_elai/glaucoma_shiley_jun26.pheno",header=T,stringsAsFactors=F,na.strings=c("NA","#N/A"))

 subs <- read.table('/mnt/sdb/genetics/ancestry_pipeline_elai/tf/theo2mil_v6B_5.fam',stringsAsFactors=F,header=F)
 names(subs)[c(1:2)] <- c("FID","IID")

 subs$order <- 1:dim(subs)[1]

 anc1 <- merge(anc,subs,by=c("FID","IID"),all.y=TRUE)
 anc2 <- merge(phenotypes,anc1,by=c("FID"),all.y=TRUE)
 anc3 <- anc2[order(anc2$order),]

#Subset to glaucoma only, black subjects
 africans <- which(anc3$bestpop_oneweek=="aam" & anc3$Diagnosis == "Glaucoma")
 
#Write out nromal subjects, for HWE violation testing 
 write.table(anc3[which(anc3$bestpop_oneweek=="aam" ),c("FID","IID")],'theoafricans_forhwe.subjects',quote=F,row.names=F)
 
 
#get a covariate sheet
 anc3_cov <- anc3[africans,]

#Load overall admixture proportions (SNPweights AIMs panel estimates) for CONTROL data
 anc <- read.table('/mnt/sdb/genetics/gradyanc/GTPC/ancestry/Adam_final_GTP_pheno.predpc_oneweek.header',stringsAsFactors=F,header=T)

 subs <- fread('pts_gtpc_aam_am-qc_v6_22.raw',data.table=F,header=F)[,1:6] #subjects 1 to 1000
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

unlist_split <- function(x, ...)
{
    toret <- unlist(strsplit(x, ...) )
    return(t(toret))
}
library(plyr)

 #Note that for the genotype data, you must be sure to code such that the reference allele is the same across datasets
for (chr in  c(c(1:22)))
{

###Load local ancestry and SNP position data

##Load positions for glaucoma data
 positions <- fread(paste("/mnt/sdb/genetics/ancestry_pipeline_elai/output/theo_linear_c2_",chr,".snpinfo.txt",sep=''), header=T,data.table=F)

#Load Local ancestry for Glaucoma data
 dat <- fread(paste("/mnt/sdb/genetics/ancestry_pipeline_elai/output/theo_linear_c2_",chr,".ps21.txt",sep=''), header=F,data.table=F)

#Make dataframes for african ancestry class
 dat_afr <- dat[africans,seq(2,dim(dat)[2],by=2)]
 names(dat_afr) <- positions$rs
 
##Load positions for control data
 pos_conta <- fread(paste("output/gtp_realquad_00_",chr,".snpinfo.txt",sep=''), header=T,data.table=F)

#Load Local ancestry for control data
 dat_conta <- fread(paste("output/gtp_realquad_00_",chr,".ps21.txt",sep=''), header=F,data.table=F)
 
#Make dataframes for african ancestry class
 dat_cont_afra <- dat_conta[,seq(2,dim(dat_conta)[2],by=2)]
 
 names(dat_cont_afra) <- pos_conta$rs
 
  
##Load positions for control data
 pos_contb <- fread(paste("output/gtp_realquad_01_",chr,".snpinfo.txt",sep=''), header=T,data.table=F)

#Load Local ancestry for control data
 dat_contb <- fread(paste("output/gtp_realquad_01_",chr,".ps21.txt",sep=''), header=F,data.table=F)

#Make dataframes for african ancestry class
 dat_cont_afrb <- dat_contb[,seq(2,dim(dat_contb)[2],by=2)]
 names(dat_cont_afrb) <- pos_contb$rs
 
 
 dat_cont_afra <- dat_cont_afra[,names(dat_cont_afra) %in% pos_contb$rs]
 dat_cont_afrb <- dat_cont_afrb[,names(dat_cont_afrb) %in% pos_conta$rs]
 
 dat_cont_afr <- rbind(dat_cont_afra,dat_cont_afrb)[africans_cont,]
 
 pos_cont <- pos_conta[pos_conta$rs %in% pos_contb$rs,]

 #check to make sure things line up
 
 table(pos_cont$rs == names(dat_cont_afr))
 
 #Load other snps to exclude (such as ambiguous markers)
ambig <- read.table('tf/ambiguous_snps.txt',header=F,stringsAsFactors=F)
names(ambig)[2] <- 'rs'


##Identify markers where calls were made for both cases/controls
 positions2 <- positions[-which(positions$rs %in% ambig$rs),]
 
 pos_overlap <- positions2[which(positions2$rs %in% pos_cont$rs) ,]
 pos_cont_overlap <- pos_cont[which(pos_cont$rs %in% pos_overlap$rs),]
 
 #check that this is correct
 print(table(pos_overlap$rs == pos_cont_overlap$rs))

#Subset case and controls to only theseoverlapping positions
 dat_afr2 <- dat_afr[,which(names(dat_afr) %in% pos_overlap$rs)]
 dat_cont_afr2 <- dat_cont_afr[,which(names(dat_cont_afr) %in% pos_overlap$rs)]


#Stack case/control data frames
 dat_afruse <- rbind(dat_afr2,dat_cont_afr2)

 
#Load genotypes, filter to select rs (in proper order) and subject

 geno <- fread(paste("/mnt/sdb/genetics/ancestry_pipeline_elai/theo2mil_v6_",chr,".raw",sep=''), header=T,data.table=F)
 names(geno)[-c(1:6)] <- sapply(names(geno)[-c(1:6)],unlist_split,split="_")[1,]
 geno2 <- geno[geno$IID %in% anc3_cov$IID,]
 geno3 <- geno2[,names(geno2) %in% names(dat_afruse)]

 
 geno_cont <- fread(paste("pts_gtpc_aam_am-qc_v6_",chr,".raw",sep=''), header=T,data.table=F)
 grep ("pts_gtpc_aam_am-qc_v6_1.raw", names(geno_cont),value=TRUE)
 
 
 names(geno_cont)[-c(1:6)] <- sapply(names(geno_cont)[-c(1:6)],unlist_split,split="_")[1,]
 geno_cont2 <- geno_cont[geno_cont$IID %in% anc3_cov_cont$IID,]
 geno_cont3 <- geno_cont2[,names(geno_cont2) %in% names(dat_afruse)]

  geno3 <-  geno3[,names(geno3) == names(geno_cont3)]
  geno_cont3 <-  geno_cont3[,names(geno_cont3) %in% names(geno3)]
  
  dat_afruse <- dat_afruse[, names(dat_afruse) %in% names(geno_cont3)]
  

#Check order of rs ids
 table(names(geno_cont3) ==  names(dat_afruse))
 table(names(geno3) ==  names(dat_afruse)) # True for all but a few,oddly
 
 table(names(geno3) ==  names(dat_afruse))

 
 
#Check order of subjects
 table(geno_cont2$IID == anc3_cov_cont$IID)
 table(geno2$IID == anc3_cov$IID)

  geno_use <- rbind(geno3,geno_cont3)
  
  geno_in <- rbind(round(dat_afruse,0),geno_use)
  
  #have to use genotypes with an additive coding!
  
  casevar <- c(rep(1,dim(dat_afr2)[1]),rep(0,dim(dat_cont_afr2)[1]))
  nsub <- length(casevar)
  #note:In the case of k=2, 1- african is european
  nsub
  
  library(metafor)
  
#Calculate MAFs in cases and controls
#Put in terms of the rare allele
mafcount <- function(x)
{
 maf1 <- sum(x,na.rm=T) / (2 * length(na.omit(x)))
 maf = maf1
 return(maf)
}

casemafs <- apply(geno3,2,mafcount)
conmafs <- apply(geno_cont3,2,mafcount)
overallmafs <- apply(geno_use,2,mafcount)

#Do regression of  case status on local ancestry 


 extractor <- function(x){ return(x$coefficients[2,1:2]) }
 glm_lancadj <- function(x,nsubs)
 {
  not_aa <- which(x[1:nsubs] == 0)
  one_aa <- which(x[1:nsubs] == 1)
  two_aa <- which(x[1:nsubs] == 2)
 


  geno_input <- x[(nsubs+1):length(x)]
  genofactor <- factor(geno_input,levels=c(0,1,2))#For construction of tables
  
  notaatable <- table(casevar[not_aa],genofactor[not_aa])
  oneaatable <- table(casevar[one_aa],genofactor[one_aa])
  twoaatable <- table(casevar[two_aa],genofactor[two_aa])

  
  notaasize <- min(notaatable )
  oneaasize <- min(oneaatable )
  twoaasize <- min(twoaatable )
  if (notaasize >= 5 )
  {
   t2_notAA  <- extractor(summary(glm(casevar[not_aa] ~ geno_input[not_aa] + ancov[not_aa,]$PC1 + ancov[not_aa,]$PC2 + ancov[not_aa,]$PC3,family='binomial' )))
  } else {t2_notAA <- NA}
  
  if (oneaasize >= 5)
  {
   t2_oneAA  <- extractor(summary(glm(casevar[one_aa] ~ geno_input[one_aa] + ancov[one_aa,]$PC1 + ancov[one_aa,]$PC2 + ancov[one_aa,]$PC3,family='binomial' )))
  } else {t2_oneAA <- NA}
  
  if (twoaasize >= 5 )
  {
   t2_twoAA  <- extractor(summary(glm(casevar[two_aa] ~ geno_input[two_aa] + ancov[two_aa,]$PC1 + ancov[two_aa,]$PC2 + ancov[two_aa,]$PC3,family='binomial' )))
  } else {t2_twoAA <- NA ; t2_oneAA <- NA }
  
  # if (notaasize >= 5 & oneaasize >= 5 & twoaasize >= 5)
  # {
   # t2_notAA  <- extractor(summary(glm(casevar[not_aa] ~ geno_input[not_aa] + ancov[not_aa,]$PC1 + ancov[not_aa,]$PC2 + ancov[not_aa,]$PC3,family='binomial' )))
   # t2_oneAA  <- extractor(summary(glm(casevar[one_aa] ~ geno_input[one_aa] + ancov[one_aa,]$PC1 + ancov[one_aa,]$PC2 + ancov[one_aa,]$PC3,family='binomial' )))
   # t2_twoAA  <- extractor(summary(glm(casevar[two_aa] ~ geno_input[two_aa] + ancov[two_aa,]$PC1 + ancov[two_aa,]$PC2 + ancov[two_aa,]$PC3,family='binomial' )))
  # } else {
   # t2_notAA <- NA
   # t2_oneAA <- NA
   # t2_twoAA <- NA
  # }
  
  
  meta_se <-  sum(t2_notAA[2]^-2 , t2_oneAA[2]^-2, t2_twoAA[2]^-2,na.rm=T)^-1
  meta_B <- sum(t2_notAA[1]*t2_notAA[2]^-2 , t2_oneAA[1]*t2_oneAA[2]^-2 , t2_twoAA[1]*t2_twoAA[2]^-2,na.rm=T) * meta_se
  
   bsa <- c(t2_notAA[1],t2_oneAA[1],t2_twoAA[1])
   sea <- c(t2_notAA[2],t2_oneAA[2],t2_twoAA[2])
    # print(table(casevar,x[1:nsubs]))
 
       # print(notaatable)
    # print(oneaatable)
    # print(twoaatable)
  # print(bsa)
   # print(sea)
   # print(meta_B)
   # print(sqrt(meta_se))
   # print(rma(yi=bsa,sei=sea,method="FE"))
  return(c(meta_B,sqrt(meta_se)))
 }
 
 
 
#Get regression p-values

   regsums <- t(apply(as.matrix(geno_in),2,glm_lancadj,nsubs=nsub))
   colnames(regsums) <- c("beta","se") 
    
   pvs <- pchisq((regsums[,1]/regsums[,2])^2,1,lower.tail=F)

   results <- data.frame(cbind(subset(pos_cont_overlap[pos_cont_overlap$rs %in% names(dat_afruse),],select=c(rs,pos,chr,maf)),regsums,pvs,casemafs,conmafs,overallmafs))

   afds <- fread('/mnt/sdb/genetics/gradyanc/hwetests/chr1_aam_casecontrol.results',header=T)
   afd_delete <- subset(afds,afd >= 0.05 | is.na(beta) | pvs < 0.3     ) #| is.na(beta)
   
   results2 <- results[!(results$rs %in% afd_delete$rs) ,]
   dim(results2)/dim(results)
   
   hwes <- fread(paste('/mnt/sdb/genetics/ancestry_pipeline_elai/theo2mil_v6_',chr,'_aam_hardy.hwe',sep=''),data.table=F)
   hwes_delete <- subset(hwes,P < 5e-3)
   
   head(results2[order(results2$pvs),])
   
   results3 <- results2[!(results2$rs %in% hwes_delete$SNP),]
   
   head(results3[order(results3$pvs),])
   
   
   
   unadj_filtered <- sort(results3$pvs)
UNADJ <- -log(unadj_filtered,10)
QQ <- -log(ppoints(length(UNADJ)),10)



GCfactor= round(median(qchisq( unadj_filtered,1,lower.tail=F),na.rm=T)/.455,3)

png(paste(chr,'_qq_noadj,.png', sep=''))
par(bty='l')

plot(c(0,max(QQ)), c(0,max(UNADJ)), xlab='Expected -log10(p)', ylab='Observed -log10(p)', col='blue', cex=1.3, cex.axis=1.2, cex.lab=1.5,pch=20)
errorbars=0
if(errorbars == 1)
{
    #code for error bars
    ranks <- c(1:length(QQ))
    CIlower <- qbeta(.025, ranks, length(QQ)-ranks +1)
    CIupper <- qbeta(.975, ranks, length(QQ)-ranks +1)
    plotCIlower <- -log(CIlower,10)
    plotCIupper <- -log(CIupper,10)
    segments(x0=QQ,x1=QQ, y0=plotCIlower,y1=plotCIupper,lwd=2,col='grey')
}

abline(0,1,col='red', lwd=2)
points(QQ, UNADJ, ,pch=20,col='blue')

legend('topleft', paste('GC Lambda =', GCfactor),  bty='n', cex=1.5, xjust=1)

dev.off()

  write.table(results3, paste("results/chr",chr,"_casecontrols.results",sep=''),row.names=F)

  }
   
#error checking - the top hit is certainly a STRAND FLIP!
  
   glm_lancadj(as.matrix(geno_in)[,which(names(geno_in)  =='rs7540868')],nsubs=nsub)
   #was rs10115595 aligned to reference??
   afds[afds$rs == 'rs7540868',]
   
   
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
