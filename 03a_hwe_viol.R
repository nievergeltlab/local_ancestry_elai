
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

#Subset to normal white subjects
 europeans <- which(anc3$bestpop_oneweek=="aam" & anc3$Diagnosis %in% c("Normal"))
 
#get a covariate sheet
 anc3_cov <- anc3[europeans,]

#Load overall admixture proportions (SNPweights AIMs panel estimates) for CONTROL data
 anc <- read.table('/mnt/sdb/genetics/gradyanc/GTPC/ancestry/Adam_final_GTP_pheno.predpc_oneweek.header',stringsAsFactors=F,header=T)

 subs <- read.table('tf/pts_gtpc_mix_am-qc_v6_9.fam',stringsAsFactors=F,header=F)
 names(subs)[c(1:2)] <- c("FID","IID")

 subs$order <- 1:dim(subs)[1]

 anc2 <- merge(anc,subs,by=c("IID"),all.y=TRUE)
 anc3_cont <- anc2[order(anc2$order),][1:1000,]
 
 #add [1:500,] for the aam only analysis

 
 
 #write.table( subset(anc3_cont,bestpop_oneweek=="aam",select=c(FID.y,IID))[1:1000,],'gtp_africans.subjects',quote=F,row.names=F)
 
 
#Subset to black subjects
 europeans_cont <- which(anc3_cont$bestpop_oneweek=="aam" ) 

 #get a covariate sheet
 anc3_cov_cont <- anc3_cont[europeans_cont ,]

#Combine case/control covariates sheets (PCs) for analysis
 ancov <- rbind(anc3_cov[,c("PC1","PC2","PC3","PC4","PC5")], anc3_cov_cont[,c("PC1","PC2","PC3","PC4","PC5")])

unlist_split <- function(x, ...)
{
    toret <- unlist(strsplit(x, ...) )
    return(t(toret))
}

#Note that for the genotype data, you must be sure to code such that the reference allele is the same across datasets

mafcount <- function(x)
{
 maf1 <- sum(x,na.rm=T) / (2 * length(na.omit(x)))
 maf = maf1
 return(maf)
}

 extractor <- function(x){ return(x$coefficients[2,1:2]) }
 glm_popdif <- function(x)
 {

  geno_input <- x
  genofactor <- factor(geno_input,levels=c(0,1,2))#For construction of tables
  
  notaatable <- table(casevar,genofactor)

  #print(notaatable)
   
  notaasize <- min(notaatable[,-1 ])

  if(notaasize > 0)
  { 
   t2_notAA  <- extractor(summary(glm(casevar[] ~ geno_input + ancov$PC1 + ancov$PC2 + ancov$PC3,family='binomial' )))
  } else {  t2_notAA <- NA }

   bsa <- c(t2_notAA[1])
   sea <- c(t2_notAA[2])

  return(c(bsa,sqrt(sea)))
 }
 
 
 
for (chr in  c(1:22))
{

 geno <- fread(paste("/mnt/sdb/genetics/ancestry_pipeline_elai/theo2mil_v6_",chr,".raw",sep=''), header=T,data.table=F)
 geno2 <- geno[geno$IID %in% anc3_cov$IID,]

 geno_cont <- fread(paste("pts_gtpc_aam_am-qc_v6_",chr,".raw",sep=''), header=T,data.table=F)
 #names(geno_cont)[-c(1:6)] <- sapply(names(geno_cont)[-c(1:6)],unlist_split,split="_")[1,]
 geno_cont2 <- geno_cont[geno_cont$IID %in% anc3_cov_cont$IID,]

 geno3 <- geno2[,names(geno2) %in% names(geno_cont2)]
 geno_cont3 <- geno_cont2[,names(geno_cont2) %in% names(geno3)]
 

 
#Check order of subjects
 table(geno_cont2$IID == anc3_cov_cont$IID)
 table(geno2$IID == anc3_cov$IID)

  geno_in <- rbind(geno3,geno_cont3)
  
  #have to use genotypes with an additive coding!
  
  casevar <- c(rep(1,dim(geno2)[1]),rep(0,dim(geno_cont2)[1]))
  nsub <- length(casevar)
  #note:In the case of k=2, 1- eurican is european
  
  #library(metafor)
  

#Do regression of  case status on local ancestry 

casemafs <- apply(geno3[,-c(1:6)],2,mafcount)
conmafs <- apply(geno_cont3[,-c(1:6)],2,mafcount)
overallmafs <- apply(geno_in[,-c(1:6)],2,mafcount)


 
# mafcount(geno_in[,which(names(geno_in) =='rs10964236_G')])
# mafcount(geno3[,which(names(geno_in) =='rs10964236_G')])
# mafcount(geno_cont3[,which(names(geno_in) =='rs10964236_G')])

 
#Get regression p-values


   regsums <- t(apply(as.matrix(geno_in)[,-c(1:6)],2,glm_popdif))
   colnames(regsums) <- c("beta","se") 
    
   pvs <- pchisq((regsums[,1]/regsums[,2])^2,1,lower.tail=F)

   results <- data.frame(regsums,pvs,casemafs,conmafs,overallmafs)
   results$afd <- abs(results$casemafs - results$conmafs)
   
   head(results[order(results$pvs),])
   results$rs <-   sapply(row.names(results),unlist_split,split="_")[1,]
   
   write.table(results, paste("hwetests/chr",chr,"_aam_casecontrol.results",sep=''),row.names=F)

   }
   
   # results[results$rs %in% badrs,]
   
   # #results2 <- merge(results,hwe,by='rs')
   # results3 <- subset(results2,P >= 5e-8 & casemafs >= .2 & conmafs > 0.2 & overallmafs > 0.2)
   
   # head(results3[order(results3$pvs),])
   
   
   # results_2 <- subset(results,casemafs > 0.05 & conmafs >= 0.05 & overallmafs >= 0.05)
  
# #error checking - the top hit is certainly a STRAND FLIP!
  
   # glm_popdif(as.matrix(geno_in)[,which(names(geno_in) =='rs10964236_G')])
   # #was rs10115595 aligned to reference??
   

# #Results
    # outdf <- data.frame(cbind(subset(pos_overlap,select=c(rs,pos,chr,maf)),regsums))
    # write.table(outdf, paste("results/chr",chr,"_caseonly_samplesite.results",sep=''),row.names=F)

# ##Plot results

# #Average ancestry per allele (Y axis)
     # avg_eurican_case <- apply(dat_eur2,2,mean,na.rm=T)
     # avg_eurican_cont <- apply(dat_cont_eur2,2,mean,na.rm=T)
     
    # # avg_european <- apply(dat_eur,2,mean,na.rm=T)

    # # mean_dif <- function(x,difval) {
        # # mean(x - difval,na.rm=T)
    # # }

    # # pos_dif <- apply(dat_eur,2,mean_dif,null_eur)
     # pos_dif <- avg_eurican_case - avg_eurican_cont

# #Note which values had significant p
     # sigvals <- which(regsums < .005)
     # sig_colors <- rep(eur_color,length(avg_eurican_case ))
    # #sig_colors[sigvals] <- "red"

# #Plot
     # pdf(paste('results/chr',chr,'_withcontglaucoma.pdf',sep=''),7,7)

     # newpos <- round(pos_overlap$pos/10000,1)
     # plot(newpos,avg_eurican_case,type="l",lwd=3,col=eur_color,ylim=c(.1,1.8),xlab="Position (KB)", ylab="Average Ancestry")
     # lines(newpos,avg_eurican_cont,type="l",lwd=3,col=eur_color)

      
     # points(newpos[sigvals],avg_eurican[sigvals],pch=16,col="red")
    # # abline(h=mean(null_eur),lwd=2,col="black",lty=2)

     # #lines(newpos,avg_european,type="l",lwd=3,col=eur_color)


     # plot(newpos,pos_dif ,type="l",lwd=3,col=eur_color,ylim=c(-.2,.2),xlab="Position (KB)", ylab="Average Ancestry Deviation")
     # points(newpos[sigvals],pos_dif[sigvals],pch=16,col="red")

     # dev.off()
     # inc=inc+1
}


###Extras





#What is the average sd per allele?
sd_eurican <- apply(dat_eur,2,sd,na.rm=T)
median(sd_eurican)
median(avg_eurican)


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
names(admixes0) <- c("eur","eur")
admixes2 <- data.frame(admixes0,subs)
admixes3 <- merge(admixes2,anc,by=c("FID","IID"),all.x=TRUE)
admixes <- admixes3[order(admixes3$order),]



t1 <- summary(lm(dat_eur[,1] ~ null_eur , data=dat_eur))
t.test(dat_eur[,1],dat_eur$null_eur,paired=T)

#Average ancestry per subject
avg_subj <- apply(dat_eur,1,mean)

mean(avg_subj)

cor.test(null_eur,avg_subj)
